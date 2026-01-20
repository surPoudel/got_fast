#!/usr/bin/env python3
"""
04_combine_chunks.py  (FAST + deterministic + MATCHES R's agrep())

- Uses Polars for fast read/split/groupby/concat
- Collapses UMIs using an agrep()-like UNANCHORED (substring) Levenshtein distance
- Precomputes match lists ONCE per BC row (big speedup) while keeping EXACT same output
- Prints timestamps + stage durations

Run:
  python -u 04_combine_chunks.py \
    --outdir /path/to/Output \
    --pcr-thresh 0.75 --lev-dist 0.2 --dupcut 2 --threads 8
"""

from __future__ import annotations

import argparse
import glob
import math
import os
import time
from datetime import datetime
from multiprocessing import Pool
from typing import Any, Dict, List, Tuple

import polars as pl


# =========================
# logging
# =========================
def log(msg: str) -> None:
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)


# =========================
# R-equivalent genotype calling
# =========================
def r_concat_call(wt: int, mut: int, pcr_thresh: float) -> str:
    denom = wt + mut
    if denom == 0:
        return "AMB"
    if (wt / denom) > pcr_thresh:
        return "WT"
    if (mut / denom) > pcr_thresh:
        return "MUT"
    return "AMB"


def r_collapse_call(sum_wt: int, sum_mut: int, pcr_ratio_thresh: float) -> str:
    denom = sum_wt + sum_mut
    if denom == 0:
        return "AMB"
    pcr_ratio = max(sum_wt, sum_mut) / denom
    if pcr_ratio > pcr_ratio_thresh:
        return "WT" if sum_wt > sum_mut else "MUT"
    return "AMB"


def r_ld_threshold(pat_len: int, ld: float) -> int:
    # R agrep(max.distance=ld): if ld < 1, it's a fraction of pattern length
    if ld < 1.0:
        return int(math.floor(pat_len * ld))
    return int(ld)


# =========================
# Polars compatibility helpers
# =========================
def _with_row_index(lf: pl.LazyFrame, name: str) -> pl.LazyFrame:
    if hasattr(lf, "with_row_index"):
        return lf.with_row_index(name)
    return lf.with_row_count(name)  # older polars


def _scan_csv(path: str, sep: str = "\t") -> pl.LazyFrame:
    if hasattr(pl, "scan_csv"):
        return pl.scan_csv(
            path,
            separator=sep,
            has_header=True,
            infer_schema_length=200,
            ignore_errors=True,
        )
    return pl.read_csv(path, separator=sep, has_header=True, infer_schema_length=200).lazy()


def _collect_compat(lf: pl.LazyFrame) -> pl.DataFrame:
    try:
        return lf.collect(engine="streaming")
    except TypeError:
        return lf.collect()


# =========================
# Utilities
# =========================
def _clean_field(x: Any) -> str:
    if x is None:
        return ""
    s = str(x).strip()
    if s == '""':
        return ""
    if len(s) >= 2 and s[0] == '"' and s[-1] == '"':
        s = s[1:-1]
    return s.strip()


def write_tsv_noquotes(path: str, header: List[str], rows: List[List[str]]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as out:
        out.write("\t".join(header) + "\n")
        for r in rows:
            out.write("\t".join("" if v == '""' else v for v in r) + "\n")


def choose_input_files(outdir: str, sample: str) -> List[str]:
    """
    Match R logic:
      - Prefer <chunk>/<sample>.summTable.txt if exists
      - Else if exactly one *.summTable.txt in chunk dir, use it
    """
    files: List[str] = []
    chunk_dirs = sorted(glob.glob(os.path.join(outdir, "[0-9][0-9][0-9][0-9]")))
    for cd in chunk_dirs:
        preferred = os.path.join(cd, f"{sample}.summTable.txt")
        if os.path.exists(preferred) and os.path.getsize(preferred) > 0:
            files.append(preferred)
            continue

        cands = sorted(glob.glob(os.path.join(cd, "*.summTable.txt")))
        cands = [p for p in cands if ".concat" not in os.path.basename(p)]
        if len(cands) == 1 and os.path.getsize(cands[0]) > 0:
            files.append(cands[0])

    return files


# =========================
# Polars: explode + aggregate pairs
# =========================
def build_agg_pairs_polars(files: List[str], pcr_thresh: float) -> pl.DataFrame:
    """
    Returns one row per (BC, UMI) with summed WT/MUT/AMB + call.in.dups + bc_first order key.
    Deterministic: we sort by (bc_first, BC, UMI) and do semicolon joining in Python.
    """
    need = ["BC", "UMI", "num.WT.in.dups", "num.MUT.in.dups", "num.amb.in.dups"]

    lfs: List[pl.LazyFrame] = []
    for f in files:
        lf = (
            _scan_csv(f)
            .select([pl.col(c).cast(pl.Utf8) for c in need])
            .with_columns([pl.col(c).fill_null("").alias(c) for c in need])
        )
        lfs.append(lf)

    lf_all = pl.concat(lfs, how="vertical")

    # Normalize BC: take first before ';' (handles "BC;BC;BC" as well as "BC")
    # Split UMI and counts into list columns and explode those; BC stays scalar and is replicated.
    lf_long = (
        lf_all.with_columns([
            pl.col("BC").str.split(";").list.first().alias("BC_one"),
            pl.col("UMI").str.split(";").alias("UMI_list"),
            pl.col("num.WT.in.dups").str.split(";").alias("WT_list"),
            pl.col("num.MUT.in.dups").str.split(";").alias("MUT_list"),
            pl.col("num.amb.in.dups").str.split(";").alias("AMB_list"),
        ])
        .select(["BC_one", "UMI_list", "WT_list", "MUT_list", "AMB_list"])
        .rename({"BC_one": "BC"})
        .explode(["UMI_list", "WT_list", "MUT_list", "AMB_list"])
        .pipe(lambda x: _with_row_index(x, "gidx"))  # first appearance ordering like R
        .with_columns([
            pl.col("UMI_list").cast(pl.Utf8).alias("UMI"),
            pl.col("WT_list").cast(pl.Int64, strict=False).fill_null(0).alias("WT"),
            pl.col("MUT_list").cast(pl.Int64, strict=False).fill_null(0).alias("MUT"),
            pl.col("AMB_list").cast(pl.Int64, strict=False).fill_null(0).alias("AMB"),
        ])
        .select(["BC", "UMI", "WT", "MUT", "AMB", "gidx"])
    )

    bc_first = lf_long.group_by("BC").agg(pl.col("gidx").min().alias("bc_first"))

    denom = pl.col("WT") + pl.col("MUT")
    call_expr = (
        pl.when(denom == 0).then(pl.lit("AMB"))
        .when((pl.col("WT").cast(pl.Float64) / denom.cast(pl.Float64)) > pcr_thresh).then(pl.lit("WT"))
        .when((pl.col("MUT").cast(pl.Float64) / denom.cast(pl.Float64)) > pcr_thresh).then(pl.lit("MUT"))
        .otherwise(pl.lit("AMB"))
    )

    pairs_lf = (
        lf_long.group_by(["BC", "UMI"])
        .agg([
            pl.sum("WT").alias("WT"),
            pl.sum("MUT").alias("MUT"),
            pl.sum("AMB").alias("AMB"),
        ])
        .join(bc_first, on="BC", how="left")
        .with_columns(call_expr.alias("call.in.dups"))
        .select(["BC", "UMI", "WT", "MUT", "AMB", "call.in.dups", "bc_first"])
        .sort(["bc_first", "BC", "UMI"])
    )

    return _collect_compat(pairs_lf)


def build_concat_rows_from_pairs(pairs: pl.DataFrame) -> List[Dict[str, str]]:
    """
    Deterministic per-BC semicolon concatenation:
      - BC order by bc_first
      - UMI sorted lexicographically within BC (already sorted in pairs)
    """
    out: List[Dict[str, str]] = []
    cur_bc: str | None = None

    umis: List[str] = []
    wts: List[str] = []
    muts: List[str] = []
    ambs: List[str] = []
    calls: List[str] = []

    for r in pairs.iter_rows(named=True):
        bc = _clean_field(r["BC"])
        umi = _clean_field(r["UMI"])
        wt = str(int(r["WT"])) if r["WT"] is not None else "0"
        mut = str(int(r["MUT"])) if r["MUT"] is not None else "0"
        amb = str(int(r["AMB"])) if r["AMB"] is not None else "0"
        call = _clean_field(r["call.in.dups"])

        if cur_bc is None:
            cur_bc = bc

        if bc != cur_bc:
            out.append({
                "BC": cur_bc,
                "UMI": ";".join(umis),
                "num.WT.in.dups": ";".join(wts),
                "num.MUT.in.dups": ";".join(muts),
                "num.amb.in.dups": ";".join(ambs),
                "call.in.dups": ";".join(calls),
            })
            cur_bc = bc
            umis, wts, muts, ambs, calls = [], [], [], [], []

        umis.append(umi)
        wts.append(wt)
        muts.append(mut)
        ambs.append(amb)
        calls.append(call)

    if cur_bc is not None:
        out.append({
            "BC": cur_bc,
            "UMI": ";".join(umis),
            "num.WT.in.dups": ";".join(wts),
            "num.MUT.in.dups": ";".join(muts),
            "num.amb.in.dups": ";".join(ambs),
            "call.in.dups": ";".join(calls),
        })

    return out


# =========================
# agrep()-like matching: UNANCHORED (substring) Levenshtein with cutoff
# =========================
def sub_levenshtein_distance(pattern: str, text: str, cutoff: int) -> int:
    """
    Minimum Levenshtein distance between pattern and ANY substring of text (unanchored),
    with early cutoff. This matches R's agrep() "approximate grep" behavior for literal strings.

    Returns cutoff+1 when it exceeds cutoff.
    """
    if cutoff < 0:
        return cutoff + 1
    if pattern == "":
        return 0

    # If text is much shorter, you need at least deletions
    if len(pattern) - len(text) > cutoff:
        return cutoff + 1

    n = len(text)

    # Unanchored: allow starting anywhere at 0 cost (distance to empty pattern is 0 at any position)
    prev = [0] * (n + 1)

    for i, pc in enumerate(pattern, start=1):
        cur = [i]
        min_row = cur[0]
        for j, tc in enumerate(text, start=1):
            ins = cur[j - 1] + 1
            dele = prev[j] + 1
            sub = prev[j - 1] + (pc != tc)

            v = ins if ins < dele else dele
            if sub < v:
                v = sub
            cur.append(v)
            if v < min_row:
                min_row = v

        if min_row > cutoff:
            return cutoff + 1
        prev = cur

    best = min(prev)  # can end anywhere => unanchored
    return best if best <= cutoff else cutoff + 1


def build_match_list_0based(umis: List[str], ld: float) -> List[List[int]]:
    """Compute directed match lists once (0-based), equivalent to repeated agrep() calls."""
    n = len(umis)
    out: List[List[int]] = [[] for _ in range(n)]
    for i, pat in enumerate(umis):
        thr = r_ld_threshold(len(pat), ld)

        if thr <= 0:
            hits = [j for j in range(n) if umis[j] == pat]
            if i not in hits:
                hits.append(i)
            out[i] = sorted(set(hits))
            continue

        hits: List[int] = []
        for j, cand in enumerate(umis):
            if sub_levenshtein_distance(pat, cand, thr) <= thr:
                hits.append(j)
        if i not in hits:
            hits.append(i)
        out[i] = sorted(set(hits))
    return out


def transitive_closure_0based(match_list: List[List[int]], active: List[bool], seed0: List[int]) -> List[int]:
    seen: set[int] = set()
    frontier = [i for i in seed0 if active[i]]
    while frontier:
        nxt: List[int] = []
        for i in frontier:
            if i in seen:
                continue
            seen.add(i)
            for j in match_list[i]:
                if active[j] and j not in seen:
                    nxt.append(j)
        frontier = nxt
    return sorted(seen)


# =========================
# FAST collapse per BC (exact output)
# =========================
def list_collapse_one_row_fast_exact(
    bc: str,
    umi_s: str,
    wt_s: str,
    mut_s: str,
    amb_s: str,
    call_s: str,
    pcr_thresh: float,
    ld: float,
    dupcut: int,
) -> Dict[str, Any]:
    bc = _clean_field(bc)

    umi_s = _clean_field(umi_s)
    wt_s = _clean_field(wt_s)
    mut_s = _clean_field(mut_s)
    amb_s = _clean_field(amb_s)
    call_s = _clean_field(call_s)

    UMIs = [x.strip().strip('"') for x in umi_s.split(";")] if umi_s else []
    WT = [int(x) for x in wt_s.split(";")] if wt_s else []
    MUT = [int(x) for x in mut_s.split(";")] if mut_s else []
    AMB = [int(x) for x in amb_s.split(";")] if amb_s else []
    CALL = [x.strip().strip('"') for x in call_s.split(";")] if call_s else []

    n = len(UMIs)
    if n == 0:
        return {
            "BC": bc,
            "UMI": "",
            "num.WT.in.dups": "",
            "num.MUT.in.dups": "",
            "num.amb.in.dups": "",
            "call.in.dups": "",
            "WT.calls": 0,
            "MUT.calls": 0,
            "amb.calls": 0,
        }

    # Pad vectors
    if len(WT) < n:
        WT += [0] * (n - len(WT))
    if len(MUT) < n:
        MUT += [0] * (n - len(MUT))
    if len(AMB) < n:
        AMB += [0] * (n - len(AMB))

    if len(CALL) != n:
        CALL = [r_concat_call(WT[i], MUT[i], pcr_thresh) for i in range(n)]
    elif len(CALL) < n:
        CALL += ["AMB"] * (n - len(CALL))

    # Precompute agrep-like neighbors ONCE (big speedup; exact result preserved via active mask)
    match_list = build_match_list_0based(UMIs, ld)
    active = [True] * n

    def compute_degrees() -> Tuple[int, int, int, List[int]]:
        deg = [0] * n
        active_count = 0
        sum_deg = 0
        max_deg = 1
        for i in range(n):
            if not active[i]:
                continue
            active_count += 1
            di = 0
            for j in match_list[i]:
                if active[j]:
                    di += 1
            deg[i] = di
            sum_deg += di
            if di > max_deg:
                max_deg = di
        return active_count, sum_deg, max_deg, deg

    active_count, sum_deg, max_deg, deg = compute_degrees()

    while sum_deg > active_count:
        # R tie-break: first max degree in current order
        to_collapse = None
        for i in range(n):
            if active[i] and deg[i] == max_deg:
                to_collapse = i
                break
        if to_collapse is None:
            break

        closure = transitive_closure_0based(match_list, active, match_list[to_collapse])

        sum_wt = sum(WT[i] for i in closure)
        sum_mut = sum(MUT[i] for i in closure)
        sum_amb = sum(AMB[i] for i in closure)

        WT[to_collapse] = sum_wt
        MUT[to_collapse] = sum_mut
        AMB[to_collapse] = sum_amb
        CALL[to_collapse] = r_collapse_call(sum_wt, sum_mut, pcr_thresh)

        for i in closure:
            if i != to_collapse:
                active[i] = False

        active_count, sum_deg, max_deg, deg = compute_degrees()

    keep_idx = [i for i in range(n) if active[i]]

    UMIs2 = [UMIs[i] for i in keep_idx]
    WT2 = [WT[i] for i in keep_idx]
    MUT2 = [MUT[i] for i in keep_idx]
    AMB2 = [AMB[i] for i in keep_idx]
    CALL2 = [CALL[i] for i in keep_idx]

    # Sort by total decreasing (stable tie-break on original order)
    totals = [WT2[i] + MUT2[i] + AMB2[i] for i in range(len(UMIs2))]
    order = sorted(range(len(UMIs2)), key=lambda i: (-totals[i], i))

    UMIs2 = [UMIs2[i] for i in order]
    WT2 = [WT2[i] for i in order]
    MUT2 = [MUT2[i] for i in order]
    AMB2 = [AMB2[i] for i in order]
    CALL2 = [CALL2[i] for i in order]

    # dupcut filter on WT+MUT (R)
    keep2 = [(WT2[i] + MUT2[i]) >= dupcut for i in range(len(UMIs2))]
    UMIs_f = [u for u, k in zip(UMIs2, keep2) if k]
    WT_f = [x for x, k in zip(WT2, keep2) if k]
    MUT_f = [x for x, k in zip(MUT2, keep2) if k]
    AMB_f = [x for x, k in zip(AMB2, keep2) if k]
    CALL_f = [c for c, k in zip(CALL2, keep2) if k]

    return {
        "BC": bc,
        "UMI": ";".join(UMIs_f),
        "num.WT.in.dups": ";".join(map(str, WT_f)),
        "num.MUT.in.dups": ";".join(map(str, MUT_f)),
        "num.amb.in.dups": ";".join(map(str, AMB_f)),
        "call.in.dups": ";".join(CALL_f),
        "WT.calls": sum(1 for c in CALL_f if c == "WT"),
        "MUT.calls": sum(1 for c in CALL_f if c == "MUT"),
        "amb.calls": sum(1 for c in CALL_f if c == "AMB"),
    }


def _mp_worker(task: Tuple[int, Dict[str, str], float, float, int]) -> Tuple[int, Dict[str, Any]]:
    idx, d, pcr, ld, dupcut = task
    out = list_collapse_one_row_fast_exact(
        d["BC"],
        d["UMI"],
        d["num.WT.in.dups"],
        d["num.MUT.in.dups"],
        d["num.amb.in.dups"],
        d["call.in.dups"],
        pcr,
        ld,
        dupcut,
    )
    return idx, out


# =========================
# Main
# =========================
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True, help="OUTROOT/Output containing 0000/,0001/...")
    ap.add_argument("--sample", default="myGoT", help="preferred chunk filename: <sample>.summTable.txt")
    ap.add_argument("--pcr-thresh", type=float, required=True)
    ap.add_argument("--lev-dist", type=float, required=True)
    ap.add_argument("--dupcut", type=int, required=True)
    ap.add_argument("--threads", type=int, required=True)
    ap.add_argument("--progress-every", type=int, default=5000, help="print progress every N BC rows during collapse")
    args = ap.parse_args()

    t0 = time.perf_counter()
    log(f"[04] start  outdir={args.outdir}  pcr={args.pcr_thresh}  ld={args.lev_dist}  dupcut={args.dupcut}  threads={args.threads}")

   # I will just overide the dupcut with value = 1 because even though we use 2 the perl script seem to use 1 and this might be more important for concat and collapse step where we will use 2 or whatever user defines
     

    outdir = os.path.abspath(args.outdir)

    # ---- discover files ----
    t = time.perf_counter()
    files = choose_input_files(outdir, args.sample)
    if not files:
        raise SystemExit(f"[04] no usable chunk summTable files found under: {outdir}")
    log(f"[04] found {len(files)} chunk files in {time.perf_counter() - t:.2f}s (first 5):")
    for f in files[:5]:
        log(f"      {f}")

    # ---- polars aggregate ----
    t = time.perf_counter()
    pairs = build_agg_pairs_polars(files, args.pcr_thresh)
    log(f"[04] polars aggregated (BC,UMI) pairs: {pairs.height:,} rows in {time.perf_counter() - t:.2f}s")

    # ---- build concat rows ----
    t = time.perf_counter()
    concat_rows_dicts = build_concat_rows_from_pairs(pairs)
    log(f"[04] built concat rows: {len(concat_rows_dicts):,} BCs in {time.perf_counter() - t:.2f}s")

    # ---- write concat ----
    t = time.perf_counter()
    concat_path = os.path.join(outdir, f"{args.sample}.summTable.concat.txt")
    write_tsv_noquotes(
        concat_path,
        ["BC", "UMI", "num.WT.in.dups", "num.MUT.in.dups", "num.amb.in.dups", "call.in.dups"],
        [
            [d["BC"], d["UMI"], d["num.WT.in.dups"], d["num.MUT.in.dups"], d["num.amb.in.dups"], d["call.in.dups"]]
            for d in concat_rows_dicts
        ],
    )
    log(f"[04] wrote concat: {concat_path} in {time.perf_counter() - t:.2f}s")

    # ---- collapse ----
    t = time.perf_counter()
    log(f"[04] starting collapse on {len(concat_rows_dicts):,} BC rows ...")

    tasks: List[Tuple[int, Dict[str, str], float, float, int]] = [
        (i, d, float(args.pcr_thresh), float(args.lev_dist), int(args.dupcut))
        for i, d in enumerate(concat_rows_dicts)
    ]

    n = len(tasks)
    results: List[Dict[str, Any] | None] = [None] * n

    if args.threads <= 1:
        for idx, d, pcr, ld, dupcut in tasks:
            results[idx] = list_collapse_one_row_fast_exact(
                d["BC"], d["UMI"], d["num.WT.in.dups"], d["num.MUT.in.dups"], d["num.amb.in.dups"], d["call.in.dups"],
                pcr, ld, dupcut
            )
            if args.progress_every and (idx + 1) % args.progress_every == 0:
                log(f"[04] collapse progress: {idx + 1:,}/{n:,}")
    else:
        done = 0
        with Pool(processes=args.threads) as pool:
            for idx, out in pool.imap_unordered(_mp_worker, tasks, chunksize=200):
                results[idx] = out
                done += 1
                if args.progress_every and done % args.progress_every == 0:
                    log(f"[04] collapse progress: {done:,}/{n:,}")

    out_rows = [r for r in results if r is not None]
    log(f"[04] collapse done in {time.perf_counter() - t:.2f}s")

    # ---- write collapsed ----
    t = time.perf_counter()
    collapsed_path = os.path.join(outdir, f"{args.sample}.summTable.concat.umi_collapsed.txt")
    write_tsv_noquotes(
        collapsed_path,
        ["BC", "UMI", "num.WT.in.dups", "num.MUT.in.dups", "num.amb.in.dups", "call.in.dups", "WT.calls", "MUT.calls", "amb.calls"],
        [
            [
                _clean_field(d["BC"]),
                _clean_field(d["UMI"]),
                _clean_field(d["num.WT.in.dups"]),
                _clean_field(d["num.MUT.in.dups"]),
                _clean_field(d["num.amb.in.dups"]),
                _clean_field(d["call.in.dups"]),
                str(int(d["WT.calls"])),
                str(int(d["MUT.calls"])),
                str(int(d["amb.calls"])),
            ]
            for d in out_rows
        ],
    )
    log(f"[04] wrote collapsed: {collapsed_path} in {time.perf_counter() - t:.2f}s")

    log(f"[04] total time: {time.perf_counter() - t0:.2f}s")


if __name__ == "__main__":
    main()


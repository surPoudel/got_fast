# Fast GoT Pipeline (LSF + Python/Polars combiner)

This repository is a **drop-in, HPC-friendly wrapper** around the original IronThrone-GoT processing workflow, with a **fast Python/Polars-based chunk combiner + UMI collapsing** step.

It is designed for workflows where:
- You have many split FASTQ chunks (hundreds to thousands)
- You want to run the **IronThrone-GoT preprocessing (Perl)** per chunk on an LSF cluster
- You want to **merge and UMI-collapse** all chunk outputs quickly using **Polars**

> Note on thresholds: the GoT manuscript describes optimizing mutation calling with a minimum duplicate read threshold of 2 and a maximum mismatch ratio of 0.2 (species-mixing study). See the manuscript for details.

---

## Repository contents

- `run_got_pipeline.sh` – main entry point (CLI)
- `00_config.sh` – defaults and environment settings
- `01_shuffle_and_split.sh` – joins R1/R2, shuffles, and splits into chunk FASTQs
- `02_submit_chunks.sh` – submits per-chunk LSF jobs
- `03_wait_for_chunks.sh` – waits for chunk completion markers
- `04_combine_chunks.py` – fast Polars combiner + R-like UMI collapsing
- `got_post_process.py` – post-processing stub (hook point)
- `IronThrone_chunk.lsf.tpl` – LSF template for per-chunk jobs
- `IronThrone_batch.lsf.tpl` – LSF template for the combine job (optional)
- `IronThrone-GoT_only_preprocess` – the Perl pipeline used for per-chunk preprocessing (passed via `--perl`)

---

## Quick start

### 1) Create a conda env (recommended)

```bash
conda create -n got_fast python=3.11 -y
conda activate got_fast
pip install polars rapidfuzz python-Levenshtein numpy
```

### 2) Run

Example (matches the command you provided):

```bash
bash run_got_pipeline.sh \
  --fastqR1 /path/to/Sample_R1_001.fastq.gz \
  --fastqR2 /path/to/Sample_R2_001.fastq.gz \
  --config /path/to/file.config \
  --whitelist /path/to/whitelist.txt \
  --sample testQR \
  --outdir /path/to/analysis \
  --queue standard \
  --dupcut 2 \
  --threads 8 \
  --max_chunks 0 \
  --perl /path/to/IronThrone-GoT_only_preprocess \
  --python-bin python \
  --python-post "$(dirname "$0")/got_post_process.py" \
  --lev-dist 0.2 \
  --pcr 0.75 \
  --bclen 16 \
  --umilen 10 \
  --keepouts 2
```

---

## Outputs

Under `Output/` (inside the run directory):

- `myGoT.summTable.concat.txt` – merged concatenation
- `myGoT.summTable.concat.umi_collapsed.txt` – per-barcode UMI-collapsed table

The same two files are copied to `--outdir`.

---

## Notes

The “slow” R collapsing step uses agrep() to build UMI neighbor sets. Two important properties of that approach:
* agrep() matching semantics are not the same as standard Levenshtein on full strings. agrep() is an approximate “grep”: it effectively allows unanchored matches (pattern can match within the candidate), which can produce a different neighbor graph than a strict full-string Levenshtein distance. Small differences in “who is a neighbor” change connected components and therefore collapsing decisions.
* The collapse algorithm is greedy and order-dependent. The R code repeatedly picks the first UMI with max degree, then expands to a transitive closure and collapses. If the graph differs slightly, the greedy choices can diverge and produce different merges downstream.

In the current Python script, collapsing uses full-string Levenshtein distances (via rapidfuzz) and builds transitive components from that graph. Even if it looks close for specific barcodes (e.g., AAACCAAAGAGCTTTG), the global behavior can differ because the neighborhood definition differs from agrep().


## Requirements
- Python: 3.9+ (3.10/3.11 recommended)
- Perl: 5.28+ (used by IronThrone-GoT_only_preprocess)
- R: not required for the fast pipeline (only needed if you still run the original R combiner)
* LSF: optional (only if using the provided *.lsf.tpl templates on an LSF cluster)

## Python packages (fast pipeline)
* These are required for 04_combine_chunks.py and got_post_process.py:
- polars (fast table reads/joins/aggregation)
- numpy
- rapidfuzz (Levenshtein distance / cdist in UMI collapsing)


### Performance

The combine step (`04_combine_chunks.py`) is designed to scale to thousands of chunks with streaming Polars parsing and Python multiprocessing for per-barcode UMI collapsing.

---

## License / attribution

The upstream GoT pipeline is available at:
- https://github.com/dan-landau/IronThrone-GoT
- https://pmc.ncbi.nlm.nih.gov/articles/PMC6782071/ (original paper)

This repository provides cluster wrappers and a fast Python combiner.

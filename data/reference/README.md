> Lightweight download: saved MAT records, per-trial CSVs, the detailed pair audit and the historical manifest are available only in the full archive. Paths to those records below describe the full dataset. Summary CSVs and figures are included here.

# Recorded reference results

Archived numerical results are separate from new runs in `results/`. The [paper-to-file guide](../../docs/PAPER_RESULTS.md) identifies manuscript evidence; [protocols](../../docs/PROTOCOLS.md) distinguish the principal and supplementary studies.

- `*_trials.csv`: individual trajectory metrics, including completion/failure fields.
- Summary CSVs, LaTeX tables, figures and controller/mapping audit CSVs.
- `*_design.mat`: full gains, sampled matrices, references and certificate diagnostics.
- `*_protocol.mat`, validation and selection records: saved inputs and evaluation cases.
- `manifest.json`: original run metadata and source hashes; the output path is normalized to `data/reference`.

MATLAB is required to inspect the saved controller arrays and regenerate plots. Per-block checkpoints and dense trajectories are not supplied; their per-trial metrics are consolidated in CSVs. This archive cannot resume a new run.

Do not use this directory as `run_all_experiments` output. Regenerate figures in a copy if you want to preserve the distributed archive exactly.

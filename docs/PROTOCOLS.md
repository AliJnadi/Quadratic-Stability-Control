> Lightweight download: saved MAT records, per-trial CSVs, the detailed pair audit and the historical manifest are available only in the full archive. Paths to those records below describe the full dataset. Summary CSVs and figures are included here.

# Experimental protocols

The MATLAB implementation uses a planar mass-spring-damper lattice with state-dependent force geometry and uncertain physical parameters. The stages below specify the synthesis and evaluation protocols. The manuscript experiments and supplementary control comparisons use separate protocols; their numerical results should not be pooled.

| Stage | Full protocol | Solver / dependency |
|---|---|---|
| mapping | All 15 masks, four derivative samples each; nonlinear closed-loop LQR sign check | lqr, symbolic model |
| sweeps | N=50/75/100/125/150/200; r=0.1/0.075/0.05/0.025; U=20/15/10/5%; all 15 masks. 29 settings, design seed 2, 1000 boundary/corner validation samples at seed 101, nominal and uncertain QSC/LQR | SeDuMi |
| diagnostic | U=5/10/15%; QSC design seeds 2–6, one nominal LQR per U; same 1000 boundary/corner cases; 36,000 trajectories | SeDuMi |
| comparison | N=100, r=0.05, U=10%; QSC/H-infinity seeds 2–6, unit-weight and weight-matched LQR once. 100 uniform cases at seed 22001 + 32 boundary/corner cases at seed 22002, three scenarios + unique nominal disturbance | SDPT3 |
| sampling | Uniform volume, centre-biased and boundary-enriched state radii; identical base directions and parameter draws; five seeds, 132 uncertain-regulation cases | SDPT3; comparison |
| saturation | Limits 0.30/0.75/2.00 N on total commanded components, same gains/cases | comparison |
| settling | First and worst final-error record in each scenario/set; tight repeats at 30 and 100 s | comparison |
| force_resolution | First and worst peak-force uniform uncertain case per gain; tight dt=0.005/0.0005 s repeats | comparison |
| radius_settling | First/worst records in each scenario at each original radius; tight 100/300 s repeats | sweeps |
| saturation_settling | First/worst final-error and most-clipped record in each scenario/set; tight 100 s repeat | comparison, saturation |
| stabilizability | All nominal + sampled physical pairs of each actuator family; raw/scaled ranks and PBH checks at 1e-10/1e-8/1e-6 | sweeps |
| bounded | Five saved comparison sample families, three force limits; maximize trace(X), initial radius 0.01, outer radius 0.05 | SDPT3; comparison |
| margin | Independent 5000-point validation at r and 1.5r, seeds 12346/12347 | SDPT3 |
| scaling | Square grids 2/3/4; N=25/50/75/100/150/200/300; independent 1000-point validation | SDPT3 |

## Reference data and new runs

The recorded full-run results are included in `data/reference/`. This archive contains saved gains, sampled matrices, input protocols, consolidated trial tables and summaries. New runs generate their own designs and outputs in `results/`; different solver versions or feasibility solutions may change the numbers.

Repeat checks use deterministic first and worst cases selected from the current run. The radius adapter uses the radius-sweep cases, and the saturation adapter also includes each scenario and validation set's most-clipped case. The selected records are saved in `*_selection.mat` files.

The sampling-law baseline is reevaluated using the comparison sample family and gain. Full mode produces 1,980 QSC trial records across the three sampling laws. Figures and LaTeX tables are generated outputs; their numbering is independent of manuscript numbering.

## Metrics

Original sweeps/diagnostic integrate with ode45 RelTol=1e-8, AbsTol=1e-10 until 100 s or a downward crossing of ||f||=1e-8, matching the corrected original event law. Position MSE uses trapezoidal time integration of the squared position-error norm. The original mixed-state, adaptive-row-averaged metric is saved separately and has no position-only m² interpretation.

Reviewer comparisons use a 30 s horizon and dt=0.002 s quadrature; RelTol=1e-7, AbsTol=1e-9 and MaxStep=0.05. Tight repeats use 1e-9/1e-11. A test escape threshold is retained from the reviewer simulation code. Settling requires terminal velocity and acceleration norms <=1e-6, final-10% peak velocity norm <=1e-6, and maximum coordinate position range <=1e-5.

Incremental requested force RMS is computed from K*e. Total commanded/applied RMS and peak components are labelled separately. Under clipping, actual performance-output effort uses applied force minus feedforward; the requested RMS column is not relabelled as applied effort.

H-infinity output scales are 0.05 for position and 1 N for incremental effort; the input disturbance is normalized by 1 N. The nominal disturbance is a smooth 0.1 N, 0.5 Hz sine/cosine burst on the first mass's x/y channels, from 2 to 12 s. Uncertain-case ratios include affine reference drift.

Aggregation first computes each gain's case mean, then averages those means across seeds. The between-seed SD is not a confidence interval. Failed trajectories are counted explicitly; per-gain case means use completed records. Aggregate figures exclude incomplete groups. Check counts before interpreting apparent improvements.

## Storage and resuming

Trial checkpoints contain metrics, not the dense ODE arrays. Every job saves its input cases and every synthesis saves K, X, Y, attempted samples, nominal/sample matrices and acceptance diagnostics. The same settings/source/environment signature is required for a resumed run. In-progress MAT writes use a temporary file followed by an atomic rename on the same filesystem.

Full sampled designs and input protocols are retained in `data/reference/`. Individual block checkpoints are omitted because their metrics are consolidated in `*_trials.csv`. The reference archive is not a resume directory for a new run.

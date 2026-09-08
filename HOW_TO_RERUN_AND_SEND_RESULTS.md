# How to regenerate corrected results and send them back

## What changed
`check_stability_samples(...)` was being called with `samples` (the N points
used for LMI/LQR synthesis) instead of `samples_test` (an independently drawn
set of N_t = 1000 held-out trajectories). Since `check_stability_samples.m`
internally caps its count to `min(ns, length(samples))`, the previously
reported MSE/SSE numbers were computed on (a subset of) the synthesis points
themselves, not on 1000 independent test trajectories as the paper states.

This has been fixed in:
- `Matlab Codes/Bulk_Test_Samples.m` (N sweep)
- `Matlab Codes/Bulk_Test_Radius.m` (r sweep)
- `Matlab Codes/Bulk_Test_Uncer.m` (uncertainty U sweep)
- `Matlab Codes/act_test.m` (actuator-count sweep)

Each now passes `samples_test` to all four `check_stability_samples` calls.
A comment block was added at the top of each file noting the fix.

CSV export was also added at each merge/aggregation step so you can hand back
plain data files instead of screenshots or transcribed table values.

## Steps to run (in MATLAB, with CVX/YALMIP + SeDuMi available as before)

From `Matlab Codes/`:

1. **N sweep** (Table: performance vs. number of samples)
   ```matlab
   Bulk_Test_Samples   % overwrites tests/N/N*.mat with corrected numbers
   cd tests/N
   make_all_plots       % regenerates trajectory_mse_vs_samples.{fig,png,eps},
                         % ss_error_vs_samples.{fig,png,eps}, time_vs_samples.{fig,png,eps}
                         % and writes merged_raw_N_sweep.csv + summary_N_sweep.csv
   cd ../..
   ```

2. **Radius sweep** (Table: performance vs. state sampling radius r)
   ```matlab
   Bulk_Test_Radius     % overwrites tests/R/N*.mat
   cd tests/R
   make_all_plots       % regenerates trajectory_mse_vs_radius.*, ss_error_vs_radius.*,
                         % time_vs_radius.* and writes summary_R_sweep.csv
   cd ../..
   ```

3. **Uncertainty sweep** (Table: performance vs. uncertainty level U)
   ```matlab
   Bulk_Test_Uncer      % overwrites tests/U/N*.mat
   cd tests/U
   make_all_plots       % regenerates trajectory_mse_vs_uncertinity.*, ss_error_vs_uncertinity.*,
                         % time_vs_uncertinity.* and writes summary_U_sweep.csv
   cd ../..
   ```

4. **Actuator-count sweep**
   ```matlab
   act_test             % overwrites test_forces/F*.mat
   cd test_forces
   plot_forces          % regenerates trajectory_MSE_vs_actuator_count.*,
                         % steady_state_error_vs_actuator_count.* and writes
                         % raw_actuator_results.csv + summary_actuator_sweep.csv
   cd ..
   ```

Each `make_all_plots` / `plot_forces` script globs `*.mat` in its own folder
(`tests/N`, `tests/R`, `tests/U`, `test_forces`), so make sure only the
freshly-overwritten `.mat` files for that sweep are present before running it
(they should overwrite in place by filename, but worth double-checking if you
ever changed the swept values in `config.m`).

## Also needed: rerun two plotting scripts (bug found in them too)

`tests/R/make_all_plots.m` and `tests/U/make_all_plots.m` had a separate bug:
their solver-time subplot plotted a flat **mean** line instead of the actual
per-point values, which hid a real trend (solver time increasing with `r`).
Both are now fixed. Please rerun:
```matlab
cd tests/R
make_all_plots        % regenerates time_vs_radius.png/.eps with the real per-r data
cd ../U
make_all_plots        % regenerates time_vs_uncertinity.png/.eps with the real per-U data
cd ../..
```
and send back `tests/R/time_vs_radius.png` and `tests/U/time_vs_uncertinity.png`
(the `.eps` versions too if you want them swapped in for the actual PDF export).
Everything else from these two folders was already correct and has been
incorporated into the manuscript.

## What to send back

Just these six CSV files (all plain text, easy for me to parse exactly):
- `tests/N/merged_raw_N_sweep.csv` and `tests/N/summary_N_sweep.csv`
- `tests/R/summary_R_sweep.csv`
- `tests/U/summary_U_sweep.csv`
- `test_forces/raw_actuator_results.csv` and `test_forces/summary_actuator_sweep.csv`

Optionally also the regenerated `.png`/`.eps` figures if you'd like me to
swap them directly into `qsc_project/img/`, but the CSVs alone are enough for
me to update `results_table.tex`, `summary_table.tex`, and the numeric claims
in `Simulation.tex`/`results.tex`/`Conclusion.tex`.

## New experiment scripts also now export CSV/PNG

`validate_lyapunov_margin.m` and `scaling_vs_dimension.m` previously only
saved `.mat` files. They now also save a `.png`/`.eps` figure and CSV
summaries, same as the four sweep scripts above:
- `validate_lyapunov_margin.m` writes `lyapunov_margin_vs_distance.png`,
  `validation_lyapunov_margin_points.csv` (every held-out point),
  `validation_lyapunov_margin_bins.csv` (binned by radial distance), and
  `validation_lyapunov_margin_summary.csv` (top-line fraction/mean/worst
  margin for each of the two held-out sets).
- `scaling_vs_dimension.m` writes `scaling_vs_dimension_summary.png`,
  `scaling_vs_dimension_summary.csv` (N\*(n) and cost per grid size), and
  `scaling_vs_dimension_detail.csv` (every tested N per grid size).

If you run these two (they're independent of the four sweeps above and
don't need to be run in any particular order), send back the CSVs and PNGs
they produce as well.

## One thing to double check before trusting the new numbers

`config.m` sets `seed = 1`, and `generate_samples(..., seed+100, 0)` is used
for `samples_test` with `mode = 0` (maximal-deviation sampling, not uniform-in-ball).
This matches the original script's intent (`samples_test = generate_samples(...,
seed+100, 0)` was already there, just unused) — I did not change the sampling
mode or seed, only which variable gets passed to `check_stability_samples`.
If `mode = 0` (fixed boundary magnitude) wasn't actually the intended test
distribution, that's a separate, independent thing to double check against
what the paper says about how test trajectories are drawn.

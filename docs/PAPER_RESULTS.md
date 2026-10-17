> Lightweight download: saved MAT records, per-trial CSVs, the detailed pair audit and the historical manifest are available only in the full archive. Paths to those records below describe the full dataset. Summary CSVs and figures are included here.

# Paper-to-file guide

This guide follows the final revised manuscript, **Robust Sampling-Based Control for Elastic Lattices with Parameter Uncertainty**. Named results are used rather than relying on figure/table numbers that may change during typesetting.

All paths below are relative to `data/reference/`.

| Manuscript result | Recorded source | New-run stage |
|---|---|---|
| Sample-count sweep, N=50–200 | `sweeps/summary.csv`, rows `sweep=N`; individual cases in `sweeps/N*/` | `sweeps` |
| Radius sweep | `sweeps/summary.csv`, rows `sweep=R`; individual cases in `sweeps/R*/` | `sweeps` |
| Single-design uncertainty sweep | `sweeps/summary.csv`, rows `sweep=U`; individual cases in `sweeps/U*/` | `sweeps` |
| Five-design uncertainty comparison | `diagnostic/aggregate_summary.csv`; seed means in `diagnostic/summary.csv` and `diagnostic/seed_summary.csv`; trial CSVs in `diagnostic/U050/`, `U100/`, `U150/` | `diagnostic` |
| Actuator masks/count comparison | `sweeps/summary.csv`, rows `sweep=F`; trial records in `sweeps/F01/` through `F15/` | `sweeps` |
| Computational environment: n=16/36/64, N=100 | `scaling/summary.csv`, rows `N=100`; times round to 2.3/26.0/340.2 s | `scaling` |
| Appendix A: sampling-law sensitivity | `sampling/aggregate_summary.csv`, `sampling/seed_summary.csv`, per-law trial files | `comparison` then `sampling` |

For principal sweeps, `summary.csv` uses `position_MSE_mean`, `final_position_error_mean` and `incremental_force_RMS_mean`. In aggregate tables the equivalent columns are `position_MSE_mean`, `final_error_mean` and `force_RMS_mean`. Use `scenario=uncertain_regulation` for uncertain results. The principal validation set is `boundary_corner`.

The five-design table uses QSC averages across five seed-level means and the sample standard deviation across those means. LQR has one nominal gain per uncertainty setting. It is not five independent LQR designs. All gains share normalized case draws. The two mean MSE reductions are about 56.0% and 89.5% at U=10% and 15% respectively; the corresponding force-RMS ratios are 4.56 and 7.26.

The actuator-count table first averages cases within each mask and then averages masks with the same active-mass count. Every active mass has two force channels. Different masks have different feasible nominal references.

`legacy_summary.csv` contains the original mixed-state adaptive-row-averaged metric. It is not the position-only time-weighted metric used in the revised numerical tables.

## Additional analyses supplied during review

The following records are available for inspecting the additional analyses. They are separate from the principal paper protocol.

| Analysis | Recorded folder |
|---|---|
| Unit-weight LQR, weight-matched LQR, QSC, sampled H-infinity | `comparison/` |
| Force clipping | `saturation/`, `saturation_settling/` |
| Long-horizon and integration-resolution checks | `settling/`, `radius_settling/`, `force_resolution/` |
| Linear stabilizability/rank audit | `stabilizability/pair_audit.csv` |
| Ellipsoidal command bounds | `bounded/controller_audit.csv` |
| Independent homogeneous-matrix margin checks | `margin/`, and `scaling/` validation columns |

The sampling and controller-comparison studies distinguish `uniform` and `stress` validation sets. The appendix contains both sets. A failure of the returned common Lyapunov matrix on a held-out Jacobian need not imply an unstable matrix. Conversely, stable simulations do not prove a continuous-domain certificate.

The saved `RUN_COUNTS.csv` reports trajectory-record counts, not statistically independent cases. Settling repeats and reused stress cases are deliberately correlated.

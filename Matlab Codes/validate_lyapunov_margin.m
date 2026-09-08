% validate_lyapunov_margin.m
%===========================================================
% INDEPENDENT HELD-OUT VALIDATION OF THE QUADRATIC STABILITY
% CERTIFICATE
%
% Motivation
% ----------
% The QSC synthesis LMI is enforced only at the N sampled
% linearization points (A_s,B_s) used during synthesis. This
% script quantifies how well that finite-sample certificate
% extends to a large, independently drawn set of held-out
% state--parameter points, by directly evaluating the Lyapunov
% derivative condition
%
%   M(z,rho) = (A(z,rho)+B(z,rho)*K)' * P + P * (A(z,rho)+B(z,rho)*K)
%
% at each held-out point, where P is the Lyapunov matrix
% recovered from the synthesized (K,P). The condition is
% satisfied at a point iff M(z,rho) is negative definite, i.e.
% iff lambda_max(M) < 0. We report:
%
%   - the fraction of held-out points satisfying lambda_max(M) < 0
%   - the mean and worst-case (max) value of lambda_max(M)
%     (the "stability margin": more negative is better)
%   - the same statistics restricted to increasingly large
%     sub-balls of the held-out region, to see whether the
%     margin degrades smoothly or breaks down sharply near the
%     boundary of the sampled region.
%
% This directly addresses the "important remaining scientific
% enhancement" noted in REVISION_NOTES.md.
%
% NOTE on solve_quadratic_stability.m convention:
%   The variable named "P" returned by solve_quadratic_stability
%   is the LMI change-of-variables matrix (paper's X = P^{-1}),
%   NOT the Lyapunov matrix itself. The true Lyapunov matrix is
%   P_lyap = inv(P_ret). This script performs that inversion
%   explicitly; double-check this matches your local copy of
%   solve_quadratic_stability.m if you have modified it.
%===========================================================

clc; clear; close all;

config;
verbose = 1;

% -----------------------------------------------------------
% 1. Build model and reference operating point (same pipeline
%    as Bulk_Test_*.m)
% -----------------------------------------------------------
model = f_grid_symbolic(row, col, force_selector);
[z_eq, report] = validate_symbolic_model(model, row, col, verbose);

[z_ref, u_ref, u_ref_full, info_op] = generate_operating_point( ...
    model, z_eq, rho_nom, u_mag, opts, eps_f, eps_z, force_selector, seed, verbose, row, col, l);

% -----------------------------------------------------------
% 2. Synthesize QSC at the nominal design point (N, r, uncert
%    taken from config.m; edit here to check a specific design)
% -----------------------------------------------------------
samples_train = generate_samples(z_ref, u_ref, rho_nom, uncert, N, r, seed+1, 1);
grid_sys      = compute_system_matrices(model, z_ref, u_ref, rho_nom, samples_train);

[n, m] = size(grid_sys.B_r);
Q = Q_eps * eye(n);
[K_QSC, P_ret, info_QSC] = solve_quadratic_stability(grid_sys, Q, verbose, 'QSC');

if ~info_QSC.feasible
    error('QSC synthesis was infeasible for the design point in config.m; cannot validate.');
end

P_lyap = inv(P_ret);
P_lyap = (P_lyap + P_lyap')/2;   % symmetrize to remove numerical asymmetry

% -----------------------------------------------------------
% 3. Held-out validation set
%    - Different RNG seed from synthesis (seed+1) and from the
%      standard test set (seed+100), so this set is independent
%      of both.
%    - Larger than N by default (M_val >> N) for a meaningful
%      independent check.
%    - Radius r_val >= r lets us probe both inside and slightly
%      beyond the synthesis region.
% -----------------------------------------------------------
M_val   = 5000;              % number of held-out points
r_val   = r;                 % state radius for held-out points (m)
r_val_ext = 1.5 * r;         % extended radius to probe beyond the synthesis ball
uncert_val = uncert;         % parameter uncertainty for held-out points

val_seed = seed + 12345;     % must not overlap with training/testing seeds

samples_val     = generate_samples(z_ref, u_ref, rho_nom, uncert_val, M_val, r_val,     val_seed,   1);
samples_val_ext = generate_samples(z_ref, u_ref, rho_nom, uncert_val, M_val, r_val_ext, val_seed+1, 1);

% -----------------------------------------------------------
% 4. Evaluate the Lyapunov derivative condition at every
%    held-out point
% -----------------------------------------------------------
[margin, radial_dist] = lyap_margin_on_samples(model, samples_val, K_QSC, P_lyap, z_ref);
[margin_ext, radial_dist_ext] = lyap_margin_on_samples(model, samples_val_ext, K_QSC, P_lyap, z_ref);

report_margin_stats('Held-out set (r_val = r)', margin);
report_margin_stats(sprintf('Held-out set (r_val = %.2f * r, extrapolation check)', r_val_ext/r), margin_ext);

% -----------------------------------------------------------
% 5. Margin vs. radial distance from the reference point, to
%    see whether degradation is smooth or abrupt near the
%    boundary of the sampled ball.
% -----------------------------------------------------------
n_bins = 10;
edges  = linspace(0, r_val_ext, n_bins+1);
bin_mean_margin = nan(n_bins,1);
bin_worst_margin = nan(n_bins,1);
bin_frac_ok = nan(n_bins,1);

for b = 1:n_bins
    idx = radial_dist_ext >= edges(b) & radial_dist_ext < edges(b+1);
    if any(idx)
        bin_mean_margin(b)  = mean(margin_ext(idx));
        bin_worst_margin(b) = max(margin_ext(idx));
        bin_frac_ok(b)      = mean(margin_ext(idx) < 0);
    end
end

figure;
subplot(2,1,1);
bin_centers = (edges(1:end-1)+edges(2:end))/2;
plot(bin_centers, bin_mean_margin, '-o', 'DisplayName', 'mean margin'); hold on;
plot(bin_centers, bin_worst_margin, '-s', 'DisplayName', 'worst-case (max) margin');
yline(0, '--k', 'stability boundary');
xline(r, ':r', 'synthesis radius r');
xlabel('Radial distance from reference state [m]');
ylabel('\lambda_{max}(M)  (< 0 required)');
legend('Location','best');
title('Lyapunov derivative margin vs. distance from reference point');
grid on;

subplot(2,1,2);
plot(bin_centers, bin_frac_ok, '-o');
yline(1, '--k');
xline(r, ':r', 'synthesis radius r');
xlabel('Radial distance from reference state [m]');
ylabel('Fraction satisfying \lambda_{max}(M) < 0');
ylim([0 1.05]);
grid on;

% -----------------------------------------------------------
% 5b. Save the figure and export plain CSV/text summaries so
%     results can be handed back without needing MATLAB to
%     re-open a .mat/.fig file.
% -----------------------------------------------------------
savefig('lyapunov_margin_vs_distance.fig');
print('-dpng', '-painters', 'lyapunov_margin_vs_distance.png');
print('-depsc2', '-painters', 'lyapunov_margin_vs_distance.eps');

% Per-point CSV for both held-out sets (r_val = r, and extended radius)
T_points = table( ...
    [repmat({'r_val=r'}, numel(margin), 1); repmat({'r_val=1.5r'}, numel(margin_ext), 1)], ...
    [radial_dist; radial_dist_ext], ...
    [margin; margin_ext], ...
    'VariableNames', {'Set','RadialDistance','LyapunovMargin'});
writetable(T_points, 'validation_lyapunov_margin_points.csv');
fprintf('Wrote validation_lyapunov_margin_points.csv (%d rows)\n', height(T_points));

% Binned summary CSV (mean/worst margin and fraction-OK vs. radial bin)
T_bins = table(bin_centers(:), bin_mean_margin(:), bin_worst_margin(:), bin_frac_ok(:), ...
    'VariableNames', {'BinCenter_m','MeanMargin','WorstCaseMargin','FractionSatisfying'});
writetable(T_bins, 'validation_lyapunov_margin_bins.csv');
fprintf('Wrote validation_lyapunov_margin_bins.csv (%d rows)\n', height(T_bins));

% Top-line scalar summary (both held-out sets), one row each
T_summary = table( ...
    {'r_val=r'; 'r_val=1.5r'}, ...
    [numel(margin); numel(margin_ext)], ...
    [mean(margin<0); mean(margin_ext<0)], ...
    [mean(margin); mean(margin_ext)], ...
    [max(margin); max(margin_ext)], ...
    [min(margin); min(margin_ext)], ...
    'VariableNames', {'Set','NumPoints','FractionSatisfying','MeanMargin','WorstCaseMargin','BestCaseMargin'});
writetable(T_summary, 'validation_lyapunov_margin_summary.csv');
fprintf('Wrote validation_lyapunov_margin_summary.csv (%d rows) — send this back along with the two CSVs above and the .png figure.\n', height(T_summary));

save('validation_lyapunov_margin_results.mat', ...
    'margin', 'radial_dist', 'margin_ext', 'radial_dist_ext', ...
    'bin_centers', 'bin_mean_margin', 'bin_worst_margin', 'bin_frac_ok', ...
    'K_QSC', 'P_lyap', 'N', 'r', 'uncert', 'M_val', 'r_val_ext');

fprintf('\nSaved: validation_lyapunov_margin_results.mat\n');

%===========================================================
% Local functions
%===========================================================
function [margin, radial_dist] = lyap_margin_on_samples(model, samples, K, P, z_ref)
    Mv = numel(samples);
    margin = zeros(Mv,1);
    radial_dist = zeros(Mv,1);

    for s = 1:Mv
        z_s   = samples(s).z;
        u_s   = samples(s).u;
        rho_s = samples(s).rho;

        A_s = model.fun.A(z_s, u_s, rho_s);
        B_s = model.fun.B_mod(z_s, u_s, rho_s);

        A_cl = A_s + B_s*K;
        Msym = A_cl' * P + P * A_cl;
        Msym = (Msym + Msym')/2;

        margin(s) = max(eig(Msym));   % negative-definite <=> max eig < 0
        radial_dist(s) = norm(z_s - z_ref);
    end
end

function report_margin_stats(label, margin)
    frac_ok = mean(margin < 0);
    fprintf('\n=== %s ===\n', label);
    fprintf('  # points                        : %d\n', numel(margin));
    fprintf('  Fraction satisfying condition   : %.4f\n', frac_ok);
    fprintf('  Mean margin  (want << 0)        : %.3e\n', mean(margin));
    fprintf('  Worst-case (max) margin         : %.3e\n', max(margin));
    fprintf('  Best-case (min) margin          : %.3e\n', min(margin));
end

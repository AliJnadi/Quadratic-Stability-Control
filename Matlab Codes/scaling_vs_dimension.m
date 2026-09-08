% scaling_vs_dimension.m
%===========================================================
% SAMPLE-COUNT SCALING VS. STATE/PARAMETER DIMENSION
%
% Motivation
% ----------
% Section 3 (Problem Statement) notes that the number of LMI
% constraints grows linearly with N, but that this does not
% imply the *required* N (to reach a given robustness margin)
% is independent of problem size. This script empirically
% measures, for lattices of increasing size (row x col), the
% smallest number of synthesis samples N*(n) that achieves a
% prescribed held-out stability margin, together with the
% associated synthesis time. This directly supports the
% "Limitations" discussion (dimensional scalability) and
% could become a new figure/table (e.g. "N* vs. state
% dimension n" and "synthesis time vs. n") if the trend is
% informative.
%
% Method
% ------
% For each lattice size in GRID_SIZES:
%   1. Build the model and reference operating point.
%   2. For N in a candidate list (coarse-to-fine or bisection):
%        a. Synthesize QSC with N samples.
%        b. Evaluate the Lyapunov margin (see
%           validate_lyapunov_margin.m) on a fixed, independent
%           held-out set of size M_VAL.
%        c. Record whether the worst-case held-out margin is
%           below the target threshold MARGIN_TARGET (more
%           negative is better; MARGIN_TARGET < 0).
%   3. N*(n) is taken as the smallest tested N for which the
%      worst-case held-out margin meets the target.
%
% This is necessarily a coarse empirical scan (SDP synthesis is
% the dominant cost), so keep GRID_SIZES and the N candidate
% list modest for a first run and refine once the trend is
% visible.
%
% WARNING: for larger lattices the symbolic model construction
% in f_grid_symbolic.m and the resulting SDP can become
% significantly slower; consider running this unattended /
% overnight for the larger grid sizes, and consider reducing
% N_CANDIDATES or M_VAL if turnaround time is a concern.
%===========================================================

clc; clear; close all;

config;   % provides rho_nom-building blocks, r, uncert, Q_eps, seed, u_mag, opts, eps_f, eps_z
verbose = 0;

% -----------------------------------------------------------
% Experiment settings (edit as needed)
% -----------------------------------------------------------
% KNOWN LIMITATION: non-square lattices (row ~= col) currently trigger
% an error elsewhere in the codebase (reported when testing this script
% with GRID_SIZES including [2 3]). Root cause not yet isolated. Until
% fixed, restrict GRID_SIZES to square lattices (row == col) only.
GRID_SIZES    = [2 2; 3 3; 4 4];   % [row col] pairs -> state dim n = 4*row*col
N_CANDIDATES  = [25 50 75 100 150 200 300];
M_VAL         = 1000;               % held-out validation points per (grid, N)
MARGIN_TARGET = -1e-3;              % require worst-case lambda_max(M) <= this
force_selector_full = @(k) ones(1,k);  % all actuators active, adjust per grid size

results = struct('row',{},'col',{},'n_state',{},'N_star',{}, ...
                  'synth_time_star',{},'worst_margin_star',{}, ...
                  'N_tested',{},'worst_margin_tested',{},'synth_time_tested',{});

for g = 1:size(GRID_SIZES,1)
    row_g = GRID_SIZES(g,1);
    col_g = GRID_SIZES(g,2);
    k_actuators = row_g*col_g;
    fs = force_selector_full(k_actuators);

    fprintf('\n================ Grid %dx%d (n_state = %d) ================\n', ...
        row_g, col_g, 4*row_g*col_g);

    m_val_local = ones(row_g*col_g,1);
    rho_nom_local = [l; k; mu_l; kd; mu_d; m_val_local];

    model = f_grid_symbolic(row_g, col_g, fs);
    [z_eq_g, ~] = validate_symbolic_model(model, row_g, col_g, verbose);

    [z_ref_g, u_ref_g, ~, ~] = generate_operating_point( ...
        model, z_eq_g, rho_nom_local, u_mag, opts, eps_f, eps_z, fs, seed, verbose, row_g, col_g, l);

    % Fixed, independent held-out set for this grid size (reused
    % across all N candidates so comparisons are apples-to-apples)
    val_seed_g = seed + 54321 + g;
    samples_val_g = generate_samples(z_ref_g, u_ref_g, rho_nom_local, uncert, M_VAL, r, val_seed_g, 1);

    N_tested = [];
    worst_margin_tested = [];
    synth_time_tested = [];

    N_star = NaN; synth_time_star = NaN; worst_margin_star = NaN;

    for N_try = N_CANDIDATES
        samples_train_g = generate_samples(z_ref_g, u_ref_g, rho_nom_local, uncert, N_try, r, seed+1, 1);
        grid_sys_g = compute_system_matrices(model, z_ref_g, u_ref_g, rho_nom_local, samples_train_g);

        [n_g, ~] = size(grid_sys_g.B_r);
        Q_g = Q_eps * eye(n_g);

        [K_g, P_ret_g, info_g] = solve_quadratic_stability(grid_sys_g, Q_g, verbose, 'QSC');

        if ~info_g.feasible
            fprintf('  N = %4d : SYNTHESIS INFEASIBLE, skipping margin check\n', N_try);
            N_tested(end+1) = N_try; %#ok<*SAGROW>
            worst_margin_tested(end+1) = NaN;
            synth_time_tested(end+1) = info_g.elapsed_time;
            continue;
        end

        P_lyap_g = inv(P_ret_g);
        P_lyap_g = (P_lyap_g + P_lyap_g')/2;

        worst_margin = -inf;
        for s = 1:numel(samples_val_g)
            z_s   = samples_val_g(s).z;
            u_s   = samples_val_g(s).u;
            rho_s = samples_val_g(s).rho;

            A_s = model.fun.A(z_s, u_s, rho_s);
            B_s = model.fun.B_mod(z_s, u_s, rho_s);
            A_cl = A_s + B_s*K_g;
            Msym = A_cl'*P_lyap_g + P_lyap_g*A_cl;
            Msym = (Msym+Msym')/2;
            worst_margin = max(worst_margin, max(eig(Msym)));
        end

        fprintf('  N = %4d : synth time = %6.3f s, worst-case held-out margin = %.3e\n', ...
            N_try, info_g.elapsed_time, worst_margin);

        N_tested(end+1) = N_try;
        worst_margin_tested(end+1) = worst_margin;
        synth_time_tested(end+1) = info_g.elapsed_time;

        if isnan(N_star) && worst_margin <= MARGIN_TARGET
            N_star = N_try;
            synth_time_star = info_g.elapsed_time;
            worst_margin_star = worst_margin;
        end
    end

    if isnan(N_star)
        fprintf('  --> Target margin %.1e NOT reached for any tested N; increase N_CANDIDATES.\n', MARGIN_TARGET);
    else
        fprintf('  --> N*(n=%d) = %d  (synth time %.3f s, worst margin %.3e)\n', ...
            4*row_g*col_g, N_star, synth_time_star, worst_margin_star);
    end

    results(end+1) = struct( ...
        'row', row_g, 'col', col_g, 'n_state', 4*row_g*col_g, ...
        'N_star', N_star, 'synth_time_star', synth_time_star, ...
        'worst_margin_star', worst_margin_star, ...
        'N_tested', N_tested, 'worst_margin_tested', worst_margin_tested, ...
        'synth_time_tested', synth_time_tested); %#ok<SAGROW>
end

% -----------------------------------------------------------
% Summary plot: N*(n) and synthesis time at N*(n) vs. n
% -----------------------------------------------------------
n_states  = [results.n_state];
N_stars   = [results.N_star];
t_stars   = [results.synth_time_star];

figure;
subplot(1,2,1);
plot(n_states, N_stars, '-o');
xlabel('State dimension n'); ylabel('N^*(n) (min. samples for target margin)');
title('Sample count needed vs. state dimension'); grid on;

subplot(1,2,2);
plot(n_states, t_stars, '-o');
xlabel('State dimension n'); ylabel('Synthesis time at N^*(n) [s]');
title('Synthesis cost at N^*(n) vs. state dimension'); grid on;

% -----------------------------------------------------------
% Save the figure and export plain CSV summaries so results can
% be handed back without needing MATLAB to re-open a .mat/.fig.
% -----------------------------------------------------------
savefig('scaling_vs_dimension_summary.fig');
print('-dpng', '-painters', 'scaling_vs_dimension_summary.png');
print('-depsc2', '-painters', 'scaling_vs_dimension_summary.eps');

% Top-line summary: one row per grid size, N*(n) and cost at N*(n)
T_summary = table( ...
    [results.row]', [results.col]', [results.n_state]', ...
    [results.N_star]', [results.synth_time_star]', [results.worst_margin_star]', ...
    'VariableNames', {'Row','Col','StateDim','N_star','SynthTime_at_Nstar_sec','WorstMargin_at_Nstar'});
writetable(T_summary, 'scaling_vs_dimension_summary.csv');
fprintf('Wrote scaling_vs_dimension_summary.csv (%d rows)\n', height(T_summary));

% Full per-(grid,N) sweep detail, one row per tested N per grid size
detail_rows = struct('Row',{},'Col',{},'StateDim',{},'N',{},'WorstMargin',{},'SynthTime_sec',{});
for gi = 1:numel(results)
    Ns_g = results(gi).N_tested;
    for ki = 1:numel(Ns_g)
        detail_rows(end+1) = struct( ... %#ok<SAGROW>
            'Row', results(gi).row, 'Col', results(gi).col, ...
            'StateDim', results(gi).n_state, 'N', Ns_g(ki), ...
            'WorstMargin', results(gi).worst_margin_tested(ki), ...
            'SynthTime_sec', results(gi).synth_time_tested(ki));
    end
end
T_detail = struct2table(detail_rows);
writetable(T_detail, 'scaling_vs_dimension_detail.csv');
fprintf('Wrote scaling_vs_dimension_detail.csv (%d rows) — send this back along with the summary CSV and the .png figure.\n', height(T_detail));

save('scaling_vs_dimension_results.mat', 'results', 'GRID_SIZES', 'N_CANDIDATES', 'M_VAL', 'MARGIN_TARGET');
fprintf('Saved: scaling_vs_dimension_results.mat\n');

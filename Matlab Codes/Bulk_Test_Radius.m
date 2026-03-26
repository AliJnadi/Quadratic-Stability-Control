clc; clear; close all;

arr = [0.1, 0.075, 0.05, 0.025];

config;
model = f_grid_symbolic(row, col, force_selector);

[z_eq, report] = validate_symbolic_model(model, row, col, verbose);

[z_ref, u_ref, u_ref_full, info] = generate_operating_point( ...
    model, z_eq, rho_nom, u_mag, opts, eps_f, eps_z, force_selector, seed, verbose, row, col, l);
verbose = 0;
for idx = arr
    r = idx;
    str = sprintf('tests/R/N%dR%dU%d.mat', N, r*1000, uncert*1000);

    samples = generate_samples(z_ref, u_ref, rho_nom, uncert, N, r, seed+1, 1);
    grid_sys = compute_system_matrices(model, z_ref, u_ref, rho_nom, samples);
    
    [n, m] = size(grid_sys.B_r);
    % Solving the problem
    Q = Q_eps * eye(n);
    [K_QSC, P_QSC, info_QSC] = solve_quadratic_stability(grid_sys, Q, verbose, 'QSC');
    [K_LQR, P_LQR, info_LQR] = solve_quadratic_stability(grid_sys, Q, verbose, 'LQR');
    
    samples_test = generate_samples(z_ref, u_ref, rho_nom, uncert, N_t, r, seed+100, 0);
    
    [info_QSC_nom, stats_QSC_nom] = check_stability_samples(model.fun.f, z_ref, u_ref, rho_nom, K_QSC, samples, [0 100], 1e-8, verbose, N_t, row, col, l, force_selector);
    [info_QSC_ran, stats_QSC_ran] = check_stability_samples(model.fun.f, z_ref, u_ref, [], K_QSC, samples, [0 100], 1e-8, verbose, N_t, row, col, l, force_selector);
    [info_LQR_nom, stats_LQR_nom] = check_stability_samples(model.fun.f, z_ref, u_ref, rho_nom, K_LQR, samples, [0 100], 1e-8, verbose, N_t, row, col, l, force_selector);
    [info_LQR_ran, stats_LQR_ran] = check_stability_samples(model.fun.f, z_ref, u_ref, [], K_LQR, samples, [0 100], 1e-8, verbose, N_t, row, col, l, force_selector);

    Final_Table;
    
    save(str, 'T')
    fprintf("Data saved %s\n", str);
end
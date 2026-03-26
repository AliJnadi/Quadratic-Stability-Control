clc; clear; close all;

config;
verbose = 0;

for i = 1:length(configs)
    bin_r = dec2bin(i, 4);
    force_selector = double(bin_r) - '0';

    str = sprintf('test_forces/F%d.mat', i);
    
    model = f_grid_symbolic(row, col, force_selector);

    [z_eq, report] = validate_symbolic_model(model, row, col, verbose);

    [z_ref, u_ref, u_ref_full, info] = generate_operating_point( ...
        model, z_eq, rho_nom, u_mag, opts, eps_f, eps_z, force_selector, seed, verbose, row, col, l); 
    
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
    
    save(str, 'T');
    fprintf("Data saved %s\n", str);
end
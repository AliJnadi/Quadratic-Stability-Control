%===========================================================
% Test Script: Nonlinear Mass–Spring Lattice System
%
% This script:
%   1) Generates the exact nonlinear grid model
%   2) Validates equilibrium correctness
%   3) Computes a reference operating point
%   4) Samples nearby states & parameters
%   5) Computes linearized system matrices for LMI design
%
% Verbose levels:
%   verbose = 0 : silent
%   verbose = 1 : pipeline progress + key norms
%   verbose = 2 : full diagnostic information
%===========================================================

clc; clear; close all;

%-----------------------------------------------------------
% Load system configuration and parameters
%-----------------------------------------------------------
config;

vprint(verbose,1,'\n=== Nonlinear Lattice System Test ===\n');

vprint(verbose,1,'Grid size               : %d x %d\n', row, col);
vprint(verbose,1,'Total masses            : %d\n', row*col);
vprint(verbose,1,'State dimension         : %d\n', 4*row*col);
vprint(verbose,1,'Sampling radius (r)     : %.3f\n', r);
vprint(verbose,1,'Number of samples (N)   : %d\n', N);
vprint(verbose,1,'Parameter uncertainty   : ±%.2f %%\n', 100*uncert);

%-----------------------------------------------------------
% Generate exact nonlinear symbolic + numeric model
%-----------------------------------------------------------
vprint(verbose,1,'\n=== Generating Exact Symbolic Model ===\n');

model = f_grid_symbolic(row, col, force_selector);

vprint(verbose,1,'Model generation completed successfully\n');

[~, m] = size(model.sym.B_mod);

vprint(verbose,1,'Number of control inputs: %d\n', m);

vprint(verbose,2,'State Jacobian size A   : %d x %d\n', ...
    size(model.sym.A,1), size(model.sym.A,2));
vprint(verbose,2,'Input Jacobian size B   : %d x %d\n', ...
    size(model.sym.B_mod,1), size(model.sym.B_mod,2));

if verbose >= 2
    vprint(verbose,2,'\nModel fields:\n');
    disp(fieldnames(model));
end

%-----------------------------------------------------------
% Validate symbolic model equilibrium structure
%-----------------------------------------------------------
vprint(verbose,1,'\n=== Validating Symbolic Model ===\n');

[z_eq, report] = validate_symbolic_model(model, row, col, verbose);

assert(report.valid, 'Symbolic grid model validation failed');

vprint(verbose,1,'Symbolic model validation PASSED\n');

%-----------------------------------------------------------
% Generate reference operating point
%-----------------------------------------------------------
vprint(verbose,1,'\n=== Generating Reference Operating Point ===\n');

[z_ref, u_ref, u_nom, info] = generate_operating_point( ...
    model, z_eq, rho_nom, row, col, delta_mag, gamma, force_selector, seed);

vprint(verbose,1,'Reference operating point generated\n');

vprint(verbose,2,'||z_ref||               : %.3e\n', norm(z_ref));
vprint(verbose,2,'||u_ref||               : %.3e\n', norm(u_ref));
vprint(verbose,1,'Residual drift norm     : %.2e\n', info.residual_norm);

%-----------------------------------------------------------
% Generate sampled states, parameters, and inputs
%-----------------------------------------------------------
vprint(verbose,1,'\n=== Generating Samples ===\n');

samples = generate_samples_with_input( ...
    model, z_ref, u_nom, rho_nom, uncert, N, r, gamma, force_selector, seed+1);

vprint(verbose,1,'%d samples (z_i, u_i, rho_i) generated\n', N);

if verbose >= 2
    vprint(verbose,2,'Sample structure fields:\n');
    disp(fieldnames(samples));
end

%-----------------------------------------------------------
% Compute linearized system matrices at reference and samples
%-----------------------------------------------------------
vprint(verbose,1,'\n=== Computing System Matrices ===\n');

grid_sys = compute_system_matrices( ...
    model, z_ref, u_ref, rho_nom, samples);

vprint(verbose,1,'Linearized system matrices computed\n');

if verbose >= 2
    vprint(verbose,2,'System matrix fields:\n');
    disp(fieldnames(grid_sys));
end

%-----------------------------------------------------------
% Finding Control Law by Lyaponuv Function and LMI
%-----------------------------------------------------------
vprint(verbose,1,'\n=== Computing Control Law ===\n');

[n, m] = size(grid_sys.B_r);
% Solving the problem
Q = Q_eps * eye(n);
[K, P, info] = solve_quadratic_stability(grid_sys, Q, verbose);

assert(info.feasible, 'Problem is not feasible');

%-----------------------------------------------------------
% Test samples (for validation)
%-----------------------------------------------------------
vprint(verbose,1,'\n=== Generating Test Samples ===\n');

samples_test = generate_samples_with_input( ...
    model, z_ref, u_nom, rho_nom, uncert, N_t, r, gamma, force_selector, seed+100);

vprint(verbose,1,'%d test samples generated\n', N_t);

results = test_closed_loop_samples(model.fun.f, K, z_ref, u_ref, P, samples_test, Tspan, force_selector, opts, verbose);

summary = show_test_results(results, z_ref, verbose);

vprint(verbose,1,'\n=== Pipeline completed successfully ===\n');

% This file is for testing non linear row by col pattern system
clc; clear; close all;

% System coiefficients
config;

vprint(verbose,1,'\n=== Nonlinear Lattice System Test ===\n');

vprint(verbose,1,'Grid size           : %d x %d\n', row, col);
vprint(verbose,1,'State dimension     : %d\n', 4*row*col);
vprint(verbose,1,'Sampling radius r   : %.3f\n', r);
vprint(verbose,1,'Number of samples N : %d\n', N);
vprint(verbose,1,'Uncertainty level   : %.2f\n', uncert);

% Generate exact symbolic model
vprint(verbose,1,'\n=== Generate Exact Symbolic Model ===\n');
model = f_grid_symbolic(row, col, force_selector);
vprint(verbose,1,'\n=== System Generated ===\n');

[~, m] = size(model.sym.B_mod);

vprint(verbose,1,'Control inputs      : %d\n', m);
vprint(verbose,2,'A matrix size       : %d x %d\n', size(model.sym.A,1), size(model.sym.A,2));
vprint(verbose,2,'B matrix size       : %d x %d\n', size(model.sym.B_mod,1), size(model.sym.B_mod,2));

if verbose >= 2
   vprint(verbose,1,'\n=== System Model ===\n\n');
   disp(fieldnames(model)) 
   vprint(verbose,1,'====================\n');
end

% f = out.f;   % exact nonlinear symbolic model
% A = out.A;   % df/dz
% B = out.B;   % df/du
% z = out.z;
% u = out.u;

% f_num = subs(f, [...], [...]);
% A0    = subs(A, [...], [...]);

% Full debugging
[z_eq, report] = validate_symbolic_model(model, row, col, verbose);

assert(report.valid, 'Symbolic grid model validation failed');

% Generating the reference point
[z_ref, u_ref, info] = generate_operating_point(model, z_eq, rho_nom, row, col, delta_mag, gamma, seed);

vprint(verbose,1,'Operating point (z_ref, u_ref) generated via least-squares drift compensation\n');
vprint(verbose,2,'Reference state norm: %.3e\n', norm(z_ref));
vprint(verbose,2,'Reference input norm: %.3e\n', norm(u_ref));
vprint(verbose,1,'Drift residual norm : %.2e\n', info.residual_norm);

% Generating sampels point (zi, u_i, \rho_i) 
vprint(verbose,1,'\n=== Generate Samples ===\n');
samples = generate_samples_with_input(model, z_ref, rho_nom, uncert, N, r, gamma, seed+1);
vprint(verbose,1,'%d Sampels (z_i, u_i, p_i) Generated', N);

if verbose >= 2
   vprint(verbose,1,'\n=== Samples Structure ===\n\n');
   disp(fieldnames(samples)) 
   vprint(verbose,1,'====================\n');
end

% Calculate A, B matrices for the reference point and sampels points
vprint(verbose,1,'\n=== Calculate System Matrices for Reference and Samples ===\n');
grid_sys = compute_system_matrices(model, z_ref, u_ref, rho_nom, samples);

if verbose >= 2
   vprint(verbose,1,'\n=== System Matrices Structure ===\n\n');
   disp(fieldnames(grid_sys)) 
   vprint(verbose,1,'===================================\n');
end
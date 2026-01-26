function [z_ref, u_ref_mod, u_nom, info] = generate_operating_point(model, z_eq, rho_nom, row, col, delta_mag, gamma, force_selector, seed)
%===========================================================
% GENERATE_OPERATING_POINT
%
% Generates a nearby operating point (z_ref, u_ref) for the
% nonlinear symbolic system:
%
%   z_dot = f(z, u, rho)
%
% using first-order Taylor expansion and regularized LS.
%
% INPUTS:
%   model     : struct from f_grid_symbolic
%   rho_nom   : nominal parameter vector (numeric or symbolic)
%   row,col   : grid size
%   delta_mag : magnitude of position perturbation
%   gamma     : LS regularization
%   seed      : RNG seed (optional)
%
% OUTPUTS:
%   z_ref     : sampled operating state
%   u_ref     : feedforward input
%   info      : diagnostic struct
%===========================================================

    % -------------------------------
    % Reproducibility
    % -------------------------------
    if nargin >= 9 && ~isempty(seed)
        rng(seed);
    end

    % -------------------------------
    % Dimensions
    % -------------------------------
    n = row;
    m = col;
    N = n * m;

    z = model.z;
    u = model.u;
    u_mod = model.u_mod;

    nu = length(u);
    nu_mod = length(u_mod);

    z_nom = double(subs(z_eq, model.params, rho_nom));
    u_nom = zeros(nu,1);
%     u_mod = zeros(nu_mod,1);

    % -------------------------------
    % Generate Reference State state (positions only)
    % -------------------------------
    delta_pos = delta_mag * (2*rand(2*N,1) - 1);
    delta = [delta_pos; zeros(2*N,1)];

    z_ref = z_nom + delta;
    
    options = optimoptions(...
                           'fsolve', ...
                           'Algorithm', 'levenberg-marquardt', ...
                           'Display', 'off');
                       
    [u_ref, ~, ~] = fsolve(@(u) model.fun.f(z_ref, u, rho_nom), u_nom, options);

    % Split results
    mask = logical([force_selector, force_selector]);
    u_nom = u_ref;
    u_ref_mod = u_nom(mask);
%     % -------------------------------
%     % Evaluate nonlinear drift f(z_ref,u_nom)
%     % -------------------------------
%     f_ref = model.fun.f(z_ref, u_nom, rho_nom);
% 
%     % -------------------------------
%     % Evaluate Jacobian B at same point
%     % -------------------------------
%     B_ref = model.fun.B_mod(z_ref, u_mod, rho_nom);
% 
%     % -------------------------------
%     % Regularized least-squares input
%     % -------------------------------
%     BtB = B_ref.' * B_ref;
%     u_ref = double(- (BtB + gamma * eye(nu_mod)) \ (B_ref.' * f_ref));

    % -------------------------------
    % Residual drift (diagnostic)
    % -------------------------------
    u_mod = u_ref;
    u_mod(mask == 0) = 0;
    residual_vec = model.fun.f(z_ref, u_mod, rho_nom);

    info.residual_norm = norm(residual_vec);
    info.delta_norm = norm(delta);
    info.seed = seed;
    info.z_nom = z_nom;
    info.u_nom = u_nom;
    info.residual_vec = residual_vec;
end
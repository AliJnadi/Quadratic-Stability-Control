function [z_ref, u_ref_mod, u_ref, info] = generate_operating_point( ...
    model, z_eq, rho_nom, u_mag, opts, eps_f, eps_z, force_selector, seed, verbose, row, col, l)
%===========================================================
% GENERATE_OPERATING_POINT
%
% Finds a nonlinear equilibrium by applying a constant random
% input and simulating the system until convergence.
%
%   z_dot = f(z, u_ref, rho)
%
% INPUTS:
%   model         : system model (with model.fun.f)
%   z_eq          : known equilibrium (usually zero-input rest)
%   rho_nom       : nominal parameters
%   u_mag         : magnitude of random reduced input
%   opts          : ODE options
%   force_selector: actuator mask
%   seed          : RNG seed (optional)
%   verbose       : 0 / 1 / 2
%
% OUTPUTS:
%   z_ref         : equilibrium state
%   u_ref_mod     : reduced input
%   u_ref         : full expanded input
%   info          : diagnostics
%===========================================================

    if nargin < 10
        verbose = 1;
    end

    if nargin >= 9 && ~isempty(seed)
        rng(seed);
    end

    % Nominal equilibrium state
    z_eq = double(subs(z_eq, model.params, rho_nom));

    % Dimensions
    nu_mod = length(model.u_mod);

    % -------------------------------
    % Random constant input
    % -------------------------------
    u_ref_mod = u_mag * (2*rand(nu_mod,1) - 1);
    u_ref     = expand_input(u_ref_mod, force_selector);

    % -------------------------------
    % Simulate nonlinear system
    % -------------------------------
    f_num = model.fun.f;

    odefun = @(t,z) f_num(z, u_ref, rho_nom);

    opts_eq = odeset(opts, 'Events', ...
        @(t,z) equilibrium_event(t, z, f_num, u_ref, rho_nom, eps_f));
    
    [t, z_traj, te, ze, ie] = ode45(odefun, [0 1000], z_eq, opts_eq);

    z_ref = z_traj(end,:)';

    % -------------------------------
    % Diagnostics
    % -------------------------------
    zdot_final = f_num(z_ref, u_ref, rho_nom);
    residual_norm = norm(zdot_final);
    
    is_equilibrium = residual_norm < eps_f;
    dz_final = z_traj(end,:)' - z_traj(end-1,:)';
    dz_residual_norm = norm(dz_final);
    state_converged = dz_residual_norm < eps_z;
    
    converged = is_equilibrium && state_converged;

    info.seed           = seed;
    info.u_ref          = u_ref;
    info.u_ref_mod      = u_ref_mod;
    info.z_traj         = z_traj;
    info.t              = t;
    info.z_ref          = z_ref;
    info.zdot_final     = zdot_final;
    info.residual_norm  = residual_norm;
    info.state_change   = dz_residual_norm;
    info.converged = converged;
    info.eps_f = eps_f;
    info.eps_z = eps_z;
    
    info.event_time = te;
    info.event_state = ze;
    info.event_index = ie;
    
    if verbose >= 1
        fprintf('  ||zdot(z_ref)|| = %.3e\n', info.residual_norm);
        fprintf('  ||Δz||          = %.3e\n', info.state_change); 
    end

    % -------------------------------
    % Optional visualization
    % -------------------------------
    if verbose >= 2
        plot_equilibrium_vs_reference(z_eq, z_ref, row, col, l);
    end
    
%     if verbose >= 2
%         % Double check the correctness of the system
%         u = zeros(size(u_ref));
%         opts_eq = odeset(opts, 'Events', ...
%             @(t,z) equilibrium_event(t, z, f_num, u, rho_nom, eps_f));
%         odefun = @(t,z) f_num(z, u, rho_nom);
%         [~, z_test, te_test, ze_test, ie_test] = ode45(odefun, [0 1000], z_ref, opts_eq);
%         
%         disp([z_eq, z_test(end, :)']);
%         disp(te_test);
%         disp(ze_test);
%         disp(ie_test);
%     end
    
end

function [value, isterminal, direction] = equilibrium_event(t, z, f_num, u_ref, rho, eps_f)
    dz = f_num(z, u_ref, rho);
    value = norm(dz) - eps_f;   % stop when <= 0
    isterminal = 1;             % stop integration
    direction = -1;             % detect decreasing
end
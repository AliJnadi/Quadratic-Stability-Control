function [info, stats] = check_stability_samples(f_num, z_ref, u_ref, rho_nom, K, samples, t_span, eps, verbose, ns,row, col, l, force_selector)
%===========================================================
% CHECK_STABILITY_SAMPLES
%
% Simulates nonlinear closed-loop system for multiple samples
% and checks whether t_span is sufficient for convergence.
%
% INPUTS:
%   f_num      : nonlinear dynamics f(z,u,rho)
%   z_ref      : reference equilibrium state
%   u_ref      : reference input
%   K          : state feedback gain
%   samples    : struct array with fields:
%                - z_s   : initial state
%                - rho_s : parameter vector
%   t_span     : [t0 tf]
%   eps        : convergence threshold for ||f||
%===========================================================

    N = length(z_ref)/4;

    ns = min(ns, length(samples));   % plot first two samples
    
    % Preallocate output
    info(ns,1) = struct( ...
        'z_traj', [], ...
        'z_des', [], ...
        'z_act', [], ...
        'error_norm', [], ...
        'dz', [], ...
        'dz_norm', [], ...
        't_final', [], ...
        'converged', []);
    
    traj_mse = zeros(ns,1);
    ss_err   = zeros(ns,1);
    
    for k = 1:ns

        z0   = samples(k).z;
        
        if isempty(rho_nom)
            rho  = samples(k).rho;
        else
            rho = rho_nom;
        end
        % -------------------------------
        % Closed-loop nonlinear dynamics
        % -------------------------------
        odefun = @(t,z) f_num( ...
            z, ...
            expand_input(u_ref + K*(z - z_ref), force_selector), ...
            rho );

        % -------------------------------
        % Event: stop when ||f|| <= eps
        % -------------------------------
        eventfun = @(t,z) convergence_event(t, z, f_num, u_ref, K, z_ref, rho, eps, force_selector);

        opts = odeset('RelTol',1e-8,'AbsTol',1e-10,'Events',eventfun);

        % -------------------------------
        % Simulate
        % -------------------------------
        [t, z_traj, te, ~, ~] = ode45(odefun, t_span, z0, opts);
        
        z_act = z_traj(end, :)';
        dz = f_num( ...
            z_act, ...
            expand_input(u_ref + K*(z_act - z_ref), force_selector), ...
            rho );
        % -------------------------------
        % Plot
        % -------------------------------
        if ns < 5
            figure;
            plot_mass_trajectories(z_traj', z_ref, row, col, l);
        end
        
        traj_mse(k) = trajectory_mse(z_traj, z_ref);
        ss_err(k)   = norm(z_ref - z_act);
        
        info(k).z_traj = z_traj;
        info(k).z_des = z_ref;
        info(k).z_act = z_act;
        info(k).error_norm = norm(z_ref - z_act);
        
        info(k).dz = dz;
        info(k).dz_norm = norm(dz);
        
        info(k).t_final = t(end);
        info(k).converged = isempty(te) == false;
    end
    
    stats.traj_mse_mean = mean(traj_mse);
    stats.traj_mse_std  = std(traj_mse);

    stats.ss_mean = mean(ss_err);
    stats.ss_std  = std(ss_err);
    
    if verbose == 2
        fprintf('\n=== Stability Diagnostic Summary ===\n');
        fprintf('%5s | %9s | %10s | %10s | %8s\n', ...
            'Sample','Converged','||z-z_ref||','||dz_final||','t_final');
        fprintf('-----------------------------------------------------\n');

        for k = 1:length(info)
            fprintf('%5d | %9d | %10.2e | %10.2e | %8.2f\n', ...
                k, info(k).converged, info(k).error_norm, info(k).dz_norm, info(k).t_final);
        end
    end

end

function [value, isterminal, direction] = convergence_event( ...
        ~, z, f_num, u_ref, K, z_ref, rho, eps, force_selector)
    
    dz = f_num(z, ...
               expand_input(u_ref + K*(z - z_ref), force_selector), ...
               rho);

    value = norm(dz) - eps;   % stop when <= 0
    isterminal = 1;           % stop integration
    direction  = -1;          % crossing from positive to zero
end

function mse = trajectory_mse(z_traj, z_ref)
    err = z_traj - z_ref.';   % broadcast
    mse = mean(sum(err.^2, 2));
end
function results = test_closed_loop_samples( ...
    f_num, K, z_ref, u_ref, P, samples_test, Tspan, force_selector, opts, verbose)
%===========================================================
% TEST_CLOSED_LOOP_SAMPLES
% Simulates closed-loop nonlinear dynamics for test samples
%
%   z_dot = f(z, u, rho)
%   u     = u_ref - K (z - z_ref)
%
% Metrics:
%   - convergence
%   - Lyapunov decrease
%   - max control effort
% 
%===========================================================

Ns = numel(samples_test);
m  = length(u_ref);

% Preallocate
results.success      = false(Ns,1);
results.trajectories = cell(Ns,1);
results.metrics(Ns,1) = struct( ...
    'final_error_norm', [], ...
    'max_error_norm', [], ...
    'max_control_norm', [], ...
    'V_decrease', [] );

for i = 1:Ns
    z0  = samples_test(i).z;
    rho = samples_test(i).rho;
    
    % Closed-loop dynamics
    odefun = @(t,z) f_num( ...
        z, ...
        expand_input(u_ref + K*(z - z_ref), force_selector), ...
        rho );

    % Simulate
    try
        [t,z] = ode45(odefun, Tspan, z0, opts);
    catch ME
        if verbose >= 1
            fprintf('ODE failed for sample %d: %s\n', i, ME.message);
        end
        continue
    end

    % Error trajectory
    e = z - z_ref.';

    % Lyapunov
    V = sum((e*P).*e,2);

    % Control effort
    u_traj = zeros(length(t), m);
    for k = 1:length(t)
        u_traj(k,:) = (u_ref + K*(z(k,:).' - z_ref)).';
    end

    % Metrics
    metrics.final_error_norm = norm(e(end,:)');
    metrics.max_error_norm   = max(vecnorm(e,2,2));
    metrics.max_control_norm = max(vecnorm(u_traj,2,2));
    metrics.V_decrease       = V(end) - V(1);

    % Success condition
    tol = 1e-3;
    results.success(i) = (metrics.final_error_norm < tol) && (metrics.V_decrease < 0);

    % Store metrics and trajectories
    results.metrics(i) = metrics;
    results.trajectories{i} = struct('t',t,'z',z,'u',u_traj,'V',V);

    if verbose >= 2
        fprintf('[Sample %d] |e(tf)|=%.2e  ΔV=%.2e\n', i, metrics.final_error_norm, metrics.V_decrease);
    end
end

results.success_rate = mean(results.success);

if verbose >= 1
    fprintf('Test success rate: %.1f %% (%d/%d)\n', ...
        100*results.success_rate, sum(results.success), Ns);
end

end

%-------------------------------------------------------
% Expand reduced input according to force_selector
%-------------------------------------------------------
function u_full = expand_input(u_reduced, force_selector)
    N = length(force_selector);       % number of masses
    u_full = zeros(2*N,1);           % full input vector

    % indices of active forces
    idx = find(force_selector);

    % x-components
    u_full(idx) = u_reduced(1:length(idx));

    % y-components
    u_full(N + idx) = u_reduced(length(idx)+1:end);
end


% function results = test_closed_loop_samples( ...
%     f_num, K, z_ref, u_ref, P, samples_test, Tspan, force_selector, opts, verbose)
% %===========================================================
% % TEST_CLOSED_LOOP_SAMPLES
% %
% % Simulates closed-loop nonlinear dynamics for test samples
% %
% %   z_dot = f(z, u, rho)
% %   u     = u_ref - K (z - z_ref)
% %
% % Metrics:
% %   - convergence
% %   - Lyapunov decrease
% %   - max control effort
% %
% % INPUTS:
% %   f_num         : function handle @(z,u,rho)
% %   K             : state feedback gain
% %   z_ref         : reference state
% %   P             : Lyapunov matrix
% %   samples_test  : struct array (z, u, rho)
% %   Tspan         : [t0 tf]
% %   opts          : ODE options
% %   verbose       : 0/1/2
% %
% % OUTPUT:
% %   results : struct with fields
% %       .success
% %       .metrics
% %       .trajectories
% %===========================================================
%     Ns = numel(samples_test);
%     m  = size(K,1);
% 
%     results.success = false(Ns,1);
%     results.metrics = struct([]);
%     results.trajectories = cell(Ns,1);
%     
%     for i = 1:Ns
%         if verbose >= 1 && mod(i,10)==0
%             fprintf('Progress: %3.0f %% (%d / %d)\n', 100*i/Ns, i, Ns);
%         end
% 
%         z0   = samples_test(i).z;
%         rho  = samples_test(i).rho;
% 
%         % Closed-loop dynamics
%         odefun = @(t,z) f_num( ...
%             z, ...
%             expand_input(u_ref, z, z_ref, K, force_selector), ...
%             rho );
% 
%         % Simulate
%         try
%             [t,z] = ode45(odefun, Tspan, z0, opts);
%         
%         catch ME
%             if verbose >= 2
%                 % Print full error immediately
%                 fprintf('❌ Sample %d failed!\n', i);
%                 fprintf('Message: %s\n', ME.message);
%                 fprintf('Identifier: %s\n', ME.identifier);
%                 fprintf('Stack trace:\n');
%             
%                 for k = 1:length(ME.stack)
%                     fprintf('  %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
%                 end
% 
%                 % Optionally stop to debug
%                 rethrow(ME);   % <-- stops execution so you can debug
%             end
%             
%             if verbose >= 1
%                 fprintf('ODE failed for sample %d\n', i);
%             end
%             continue
%         end
% 
%         % Error trajectory
%         e = z - z_ref.';
% 
%         % Lyapunov
%         V = sum((e*P).*e,2);
% 
%         % Control effort
%         u_traj = zeros(length(t), m);
%         for k = 1:length(t)
%             u_traj(k,:) = (u_ref - K*(z(k,:).' - z_ref)).';
%         end
% 
%         % Metrics
%         metrics.final_error_norm = norm(e(end,:)');
%         metrics.max_error_norm   = max(vecnorm(e,2,2));
%         metrics.max_control_norm = max(vecnorm(u_traj,2,2));
%         metrics.V_decrease       = V(end) - V(1);
% 
%         % Success condition
%         tol = 1e-3;
%         success = (metrics.final_error_norm < tol) && ...
%                   (metrics.V_decrease < 0);
% 
%         results.success(i) = success;
%         results.metrics(i) = metrics;
%         results.trajectories{i} = struct('t',t,'z',z,'u',u_traj,'V',V);
%         
%         if verbose >= 2
%             fprintf('[Sample %d] |e(tf)|=%.2e  ΔV=%.2e\n', ...
%                 i, metrics.final_error_norm, metrics.V_decrease);
%         end
%     end
% 
%     results.success_rate = mean(results.success);
% 
%     if verbose >= 1
%         fprintf('Test success rate: %.1f %% (%d/%d)\n', ...
%             100*results.success_rate, sum(results.success), Ns);
%     end
% end
% 
% %-------------------------------------------------------
% % Force selector expansion (inside ODE function)
% %-------------------------------------------------------
% function u_full = expand_input(u_reduced, z, z_ref, K, force_selector)
%     % Number of masses
%     N = length(force_selector);
%     
%     % Full input vector [ux1..uxN, uy1..uyN]
%     u_full = zeros(2*N,1);
%     
%     % Indices where force is applied
%     idx = find(force_selector);
%     
%     % Expand x-components (ux)
%     u_full(idx) = u_reduced(1:length(idx));
%     
%     % Expand y-components (uy)
%     u_full(N + idx) = u_reduced(length(idx)+1:end);
%     
%     % Add feedback if needed
%     if ~isempty(K)
%         u_full = u_full - expand_input(K*(z-z_ref), z, z_ref, [], force_selector);
%     end
% end
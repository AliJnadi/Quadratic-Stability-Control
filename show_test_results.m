function summary = show_test_results(results, z_ref, verbose)
%===========================================================
% SHOW_TEST_RESULTS
%
% Displays tables and plots for closed-loop test simulations
%
% INPUTS:
%   results : struct from test_closed_loop_samples
%   z_ref   : reference state
%   verbose : 0 = silent
%             1 = summary only
%             2 = summary + plots
%
% OUTPUT:
%   summary : table of per-sample metrics
%===========================================================

    if nargin < 3
        verbose = 1;
    end

    Ns = numel(results.metrics);

    %-----------------------------------------------------------
    % Convert metrics to table
    %-----------------------------------------------------------
    summary = struct2table(results.metrics);

    summary.Sample = (1:Ns).';
    summary.Success = results.success(:);

    summary = movevars(summary, {'Sample','Success'}, 'Before', 1);

    %-----------------------------------------------------------
    % Console summary
    %-----------------------------------------------------------
    if verbose >= 1
        fprintf('\n=== Closed-Loop Test Summary ===\n');
        fprintf('Samples tested      : %d\n', Ns);
        fprintf('Success rate        : %.1f %%\n', 100*results.success_rate);
        fprintf('Worst final error   : %.2e\n', max(summary.final_error_norm));
        fprintf('Worst max error     : %.2e\n', max(summary.max_error_norm));
        fprintf('Worst control norm  : %.2e\n', max(summary.max_control_norm));
    end

    %-----------------------------------------------------------
    % Display table (MATLAB UI)
    %-----------------------------------------------------------
    if verbose >= 1
        disp(summary);
    end

    %-----------------------------------------------------------
    % Plot worst-case trajectories
    %-----------------------------------------------------------
    if verbose >= 2
        % Worst-case based on max error
        [~, idx] = max(summary.max_error_norm);
        traj = results.trajectories{idx};

        % Error norm plot
        figure;
        plot(traj.t, vecnorm(traj.z - z_ref.',2,2), 'LineWidth', 2);
        xlabel('Time');
        ylabel('||z - z_{ref}||');
        title(sprintf('Worst-Case Error Trajectory (Sample %d)', idx));
        grid on;
        
        figure;
        plot_mass_trajectories(traj.z', z_ref, 2, 2)
        
        % Lyapunov plot
        figure;
        plot(traj.t, traj.V, 'LineWidth', 2);
        xlabel('Time');
        ylabel('V(e)');
        title(sprintf('Worst-Case Lyapunov Function (Sample %d)', idx));
        grid on;

        % Control effort plot
        figure;
        plot(traj.t, vecnorm(traj.u,2,2), 'LineWidth', 2);
        xlabel('Time');
        ylabel('||u||');
        title(sprintf('Worst-Case Control Effort (Sample %d)', idx));
        grid on;
    end

end
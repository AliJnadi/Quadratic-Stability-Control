function plot_mass_trajectories(z, z_ref, row, col)

    if ndims(z) == 3
        z = squeeze(z);
    end

    N  = row * col;
    Nt = size(z, 2);

    X  = z(1:N, :);
    Y  = z(N+1:2*N, :);

    Xr = z_ref(1:N);
    Yr = z_ref(N+1:2*N);
    
    hold on;

    % --- store handles for legend ---
    h_traj  = gobjects(1,1);
    h_start = gobjects(1,1);
    h_end   = gobjects(1,1);

    for i = 1:N
        h_traj  = plot(X(i,:), Y(i,:), 'k-', 'LineWidth', 1.2);
        h_start = plot(X(i,1), Y(i,1), 'ro', 'MarkerSize', 6);
        h_end   = plot(X(i,end), Y(i,end), 'ks', 'MarkerSize', 6);
    end

    % --- reference (single handle) ---
    h_ref = plot(Xr, Yr, 'bo', ...
        'MarkerSize', 7, ...
        'LineWidth', 2, ...
        'MarkerFaceColor', 'b');

    % --- formatting ---
    axis equal;
    grid on;
    xlabel('x');
    ylabel('y');
    title('Mass Grid Trajectories vs Reference Configuration');

    % --- explicit legend ---
    legend([h_traj, h_start, h_end, h_ref], ...
           {'Trajectory','Start','End','Reference'}, ...
           'Location','bestoutside');

end
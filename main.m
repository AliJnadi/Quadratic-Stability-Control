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

% Generate force selector and nominal model
force_selector = create_force_selector(row, col, mode);
[A_nom, B_nom, b_nom] = gen_sys_matrices(row, col, force_selector, rho_nom);
[n, m] = size(B_nom);

vprint(verbose,1,'Control inputs      : %d\n', m);
vprint(verbose,2,'A matrix size       : %d x %d\n', size(A_nom,1), size(A_nom,2));
vprint(verbose,2,'B matrix size       : %d x %d\n', size(B_nom,1), size(B_nom,2));

% Generating the reference control
[x_ref, u_ref, info] = generate_operating_point(A_nom, B_nom, b_nom, rho_nom, row, col, delta_mag, gamma, seed);
fprintf('Drift residual norm: %.2e\n', info.residual);

vprint(verbose,1,'Operating point (x_ref, u_ref) generated via least-squares drift compensation\n');
vprint(verbose,2,'Reference state norm: %.3e\n', norm(x_ref));
vprint(verbose,2,'Reference input norm: %.3e\n', norm(u_ref));
vprint(verbose,1,'Drift residual norm : %.2e\n', info.residual);

% We need to find linearisation error, by find jacobian for random points
samples = generate_samples(x_ref, rho_nom, uncert, N, r, seed);
models  = evaluate_models(samples, row, col, force_selector, model_mode);

vprint(verbose,1,'Sampled models      : %d\n', numel(models));

if verbose >= 2
    DeltaA_norms = arrayfun(@(m) norm(m.DeltaA,'fro'), models);
    vprint(verbose,2,'Jacobian Fro norm   : min %.2e / max %.2e\n', ...
        min(DeltaA_norms), max(DeltaA_norms));
end

vprint(verbose,1,'=== Setup complete ===\n\n');

vprint(verbose,1,'Model evaluation mode: %s\n', model_mode);
if strcmp(model_mode,'nominal_plus_nonlinear')
    vprint(verbose,2,'Using nominal (A,B) with sampled nonlinear Jacobians\n');
else
    vprint(verbose,2,'Using fully sampled linearized models\n');
end

% Solving the problem
Q = Q_eps * eye(n);

switch model_mode
    case 'nominal_plus_nonlinear'
        A_samples = cat(3, models.DeltaA);
        solve_quadratic_stability(A_nom, B_nom, A_samples, Q, verbose, 'nominal_plus_nonlinear');

    case 'full_sampled_linear'
        A_samples = cat(3, models.A) + cat(3, models.DeltaA);
        B_samples = cat(3, models.B); % optional
        solve_quadratic_stability(A_nom, B_nom, A_samples, B_samples, Q, verbose, 'full_sampled_linear');   
end

return

% Find control law for nonlinear system
I_n = eye(n);
Q = 0.0001*I_n;
cvx_clear;
cvx_begin SDP
    variable P(n, n) symmetric
    variable L(m, n)
    
    minimize(trace(P)*0);

    subject to:
        P >= 0;
        A*P + P*A' - B*L - L'*B' + Q <= 0;
        for i = 1:100
            (A + A_sampels(:, :, i))*P + P*(A + A_sampels(:, :, i))' - B*L - L'*B' + Q <= 0;
        end
cvx_end

if strcmp(cvx_status, 'Solved')
    disp("Found solution for nonlinear \dot x = Ax + Bu + f(x)");
    
    K = L*pinv(P);
    
    % Changing initial conditions          
    delta_x = 0.04 * (-1 + 2*rand(row*col,1));
    delta_y = 0.04 * (-1 + 2*rand(row*col,1));
    delta_dx = 0.01 * (-1 + 2*rand(row*col,1));
    delta_dy = 0.01 * (-1 + 2*rand(row*col,1));
    delta = [delta_x; delta_y; delta_dx; delta_dy];
    x0 = x_ref + delta;
    
    
    [t, x] = ode45(@(t, x) MSD_non_linear_system(t, x, x_ref, u_ref, K, A, B, b, row, col, m, muxy, kxy, lx, ly, lxy), tspan, x0);
    
    x = x';
    t = t';

    motion_matrix(:, :, 1) = x(1:row*col, :);
    motion_matrix(:, :, 2) = x(row*col + 1: 2*row*col, :);
    [r, c] = size(motion_matrix, [1 2]); 
    reference_motion(:, :, 1) = repmat(x_ref(1:row*col, :), 1, size(motion_matrix, 2));
    reference_motion(:, :, 2) = repmat(x_ref(row*col + 1: 2*row*col, :), 1, size(motion_matrix, 2));
    
    animate_mass_grid_ref_error(motion_matrix, row, col, 'mass_animation.mp4', false, reference_motion, mode)
    
    mechanical_energy = energy_calc(x, row, col, m, k, lx, ly, kxy, lxy);
    
    force = zeros(m, length(t)); 
    for i = 1:length(t) 
        [~, force(:, i)] = MSD_non_linear_system(t(i), x(:, i), x_ref, u_ref, K, A, B, b, row, col, m, muxy, kxy, lx, ly, lxy);
    end
    force_power = force_power_calc(force, x);
    damper_power = damper_power_calc(x, row, col, mu , muxy);
    
%     animate_mass_grid_ref_error_energy_power_optimized(motion_matrix, row, col, 'nl_p_Fx.mp4', true, reference_motion, mode, mechanical_energy, t, force_power, damper_power);
    
    e = x - x_ref;
    
    N = length(e);
    V = zeros(N, 1);
    dV = zeros(N, 1);
    for i = 1:N
        V(i) = e(:, i)' * P * e(:, i); 
        dV(i) = e(:, i)' * (A*P + P*A' - B*L - L'*B' + Q) * e(:, i);
    end
    
%     plotLyapunovResults(t, e, V, dV, false, false)
    return
    figure;

    % Left y-axis for V
    yyaxis left
    plot(t, V, 'b', 'LineWidth', 1.2);  % thinner line
    ylabel('V','Interpreter','latex','FontSize',16);

    % Right y-axis for dV
    yyaxis right
    plot(t, dV, 'r', 'LineWidth', 1.2);  % thinner line
    ylabel('$\dot{V}$','Interpreter','latex','FontSize',16);
    
    % Add horizontal reference line at V = 0
    yyaxis left
    yline(0,'--','Color',[0.2 0.2 0.2],'LineWidth',1);
    
    % Add horizontal reference line at dV = 0
    yyaxis right
    yline(0,'--','Color',[0.2 0.2 0.2],'LineWidth',1);

    % Common x-axis
    xlabel('Time','Interpreter','latex','FontSize',16);
    title('Lyapunov Function and its Derivative','Interpreter','latex','FontSize',18);

    % Make axes ticks larger
    ax = gca;
    ax.XAxis.FontSize = 14;
    ax.YAxis(1).FontSize = 14;  % left axis
    ax.YAxis(2).FontSize = 14;  % right axis

    % Legend
    legend('V', '$V_{ref}$', '$\dot{V}$', '$\dot{V}_{ref}$','Interpreter','latex','FontSize',14,'Location','northeast');

    % Activate grid
    grid on;
    box on;

    
    return
    
    non_linear.t = t;
    non_linear.x = x;
    
    non_linear.K = K;
    non_linear.P = P;
    non_linear.S = S;
    
    non_linear.x0 = x0;
    non_linear.x_ref = x_ref;
    non_linear.u_ref = u_ref;
    
    non_linear.P = P;
    non_linear.S = S;
    
    non_linear.motion_matrix = motion_matrix;
    non_linear.reference_motion = reference_motion;
    
    non_linear.force = force;
    non_linear.mechanical_energy = mechanical_energy;
    non_linear.force_power = force_power;
    non_linear.damper_power = damper_power;
    
    save('results_point_F_x.mat', 'non_linear');
else
    disp("No solution found."); 
end

function A_samples = estimate_F(f_handle, x_ref, pho_samples, N, r)
    m = 5;                % number of trapezoid nodes
    n = length(x_ref);

    % ----- 1. Sample error vectors inside ball -----
    D = randn(n,N);
    D = D ./ vecnorm(D);       % normalize directions
    radii = r * rand(1,N).^(1/n);
    E = D .* radii;            % sampled e vectors

    % Trapezoid integration points (midpoints inside (0,1))
    s = linspace(0,1,m+2);
    s = s(2:end-1);            % remove endpoints 0 and 1
    w = ones(size(s)) / length(s);   % uniform weights

    % ----- 2. Compute ΔF(e_i) for all samples -----
    A_samples = zeros(n,n,N);

    for i = 1:N
        e_i = E(:,i);
        DeltaF_i = zeros(n,n);

        % numerical integration: average of Jacobians
        for k = 1:length(s)
            xk = x_ref + s(k)*e_i;
            Jk = numericJacobian(f_handle, xk);
            DeltaF_i = DeltaF_i + w(k)*Jk;
        end

        % matrix for LMI
        A_samples(:,:,i) = DeltaF_i;
    end
end

function [dx, u] = MSD_non_linear_system(t, x, x_ref, u_ref, K, A, B, b, row, col, m, muxy, kxy, lx, ly, lxy)
    u = u_ref - K*(x - x_ref);
    dx = A*x + B*u + b + non_linear_sys(x, row, col, m, muxy, kxy, lx, ly, lxy);
end
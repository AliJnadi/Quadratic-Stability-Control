function models = evaluate_models(samples, row, col, force_selector, model_mode)
%EVALUATE_MODELS
% Generates linearized models according to selected approximation mode

    N = numel(samples);

    switch model_mode

        case 'nominal_plus_nonlinear'
            models(N,1) = struct('DeltaA', []);

            for s = 1:N
                x   = samples(s).x;
                rho = samples(s).rho;

                DeltaA = numericJacobian( ...
                    @(z) non_linear_sys(z, row, col, rho), x );

                models(s).DeltaA = DeltaA;
            end

        case 'full_sampled_linear'
            models(N,1) = struct('A', [], 'B', [], 'DeltaA', []);

            for s = 1:N
                x   = samples(s).x;
                rho = samples(s).rho;

                [A, B, ~] = gen_sys_matrices(row, col, force_selector, rho);

                DeltaA = numericJacobian( ...
                    @(z) non_linear_sys(z, row, col, rho), x );

                models(s).A = A;
                models(s).B = B;
                models(s).DeltaA = DeltaA;
            end

        otherwise
            error('Unknown model_mode: %s', model_mode);
    end
end

% function models = evaluate_models(samples, row, col, force_selector)
%     N = length(samples);
%     models(N) = struct('A', [], 'B', [], 'b', [], 'DeltaA', []);
% 
%     for s = 1:N
%         x = samples(s).x;
%         rho = samples(s).rho;
% 
%         % Linear part
%         [A, B, b] = gen_sys_matrices(row, col, force_selector, rho);
% 
%         % Nonlinear Jacobian
%         f_handle = @(x) non_linear_sys(x, row, col, rho);
%         DeltaA = numericJacobian(f_handle, x);
% 
%         % Store
%         models(s).A = A;
%         models(s).B = B;
%         models(s).b = b;
%         models(s).DeltaA = DeltaA;
%     end
% end
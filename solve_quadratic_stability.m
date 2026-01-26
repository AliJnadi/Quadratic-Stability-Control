function [K, P, info] = solve_quadratic_stability(system_struct, Q, verbose)
%SOLVE_QUADRATIC_STABILITY
% Solves sampled quadratic stability LMIs for a nonlinear system
%
% Linearized sampled dynamics:
%   dot(e) = A_i e + B_i u
% State feedback:
%   u = K e
%
% A common quadratic Lyapunov function V = e' P e is sought such that:
%   (A_i + B_i K)' P + P (A_i + B_i K) + Q <= 0
% for the reference model and all sampled models.
%
% INPUTS:
%   system_struct : struct with fields
%       A_r  (n x n)       reference A matrix
%       B_r  (n x m)       reference B matrix
%       A_s  (n x n x Ns)  sampled A matrices
%       B_s  (n x m x Ns)  sampled B matrices
%
%   Q        : positive definite weighting matrix (n x n)
%   verbose  : 0 = silent, 1 = normal, 2 = detailed
%
% OUTPUTS:
%   K    : state feedback gain
%   P    : Lyapunov matrix
%   info : diagnostic information struct
%
% NOTES:
% - Reference model is NOT included in samples
% - No affine term is used (already compensated beforehand)
% - Uses L = K P change of variables for convexity

    % -----------------------------
    % Defaults
    % -----------------------------
    if nargin < 3 || isempty(verbose)
        verbose = 1;
    end

    % -----------------------------
    % Extract system data
    % -----------------------------
    A_r = system_struct.A_r;
    B_r = system_struct.B_r;

    A_s = system_struct.A_s;
    B_s = system_struct.B_s;

    n  = size(A_r,1);
    m  = size(B_r,2);
    Ns = size(A_s,3);

    % -----------------------------
    % Basic checks
    % -----------------------------
    assert(size(Q,1) == n && size(Q,2) == n, 'Q must be n x n');
    assert(issymmetric(Q), 'Q must be symmetric');
    assert(size(B_s,3) == Ns, 'A_s and B_s must have same number of samples');

    % -----------------------------
    % Info struct
    % -----------------------------
    info = struct();
    info.Ns = Ns;

    % -----------------------------
    % Solve LMI
    % -----------------------------
    tic;
    cvx_clear;

    if verbose == 0
        cvx_begin sdp quiet
    else
        cvx_begin sdp
    end

        variable P(n,n) symmetric
        variable L(m,n)

        minimize(0)

        subject to
            % Lyapunov matrix
            P >= 1e-6 * eye(n);

            % Reference model constraint
            A_r*P + P*A_r' + B_r*L + L'*B_r' + Q <= 0;

            % Sampled models
            for i = 1:Ns
                A_i = A_s(:,:,i);
                B_i = B_s(:,:,i);

                A_i*P + P*A_i' + B_i*L + L'*B_i' + Q <= 0;
            end

    cvx_end
    info.elapsed_time = toc;
    info.cvx_status   = cvx_status;
    info.cvx_optval   = cvx_optval;

    % -----------------------------
    % Post-processing
    % -----------------------------
    if strcmp(cvx_status,'Solved') || strcmp(cvx_status,'Inaccurate/Solved')
        K = L / P;

        info.feasible   = true;
        info.min_eig_P  = min(eig(P));
        info.trace_P    = trace(P);
        info.norm_K     = norm(K,'fro');

        if verbose >= 1
            fprintf('✔ Quadratic stability achieved with %d samples.\n', Ns);
            fprintf('  Solve time   : %.3f s\n', info.elapsed_time);
        end

        if verbose >= 2
            fprintf('  min eig(P)   : %.3e\n', info.min_eig_P);
            fprintf('  trace(P)     : %.3e\n', info.trace_P);
            fprintf('  ||K||_F      : %.3e\n', info.norm_K);
        end
    else
        K = [];
        P = [];

        info.feasible = false;

        if verbose >= 1
            fprintf('✘ Quadratic stability NOT feasible (%s)\n', cvx_status);
        end
    end
end

function sys = compute_system_matrices(model, z_ref, u_ref, rho_nom, samples)
%===========================================================
% COMPUTE_SYSTEM_MATRICES
%
% Computes linearized system matrices at:
%   - Reference point (z_ref, u_ref, rho_nom)
%   - Sampled points (z_s, u_s, rho_s)
%
% Linearization:
%   z_dot ≈ A z + B u + c
%
% INPUTS:
%   model   : model struct with numeric handles
%             model.fun.A, model.fun.B
%   z_ref   : reference state (numeric)
%   u_ref   : reference input (numeric)
%   rho_nom : nominal parameter vector (numeric)
%   samples : struct array with fields
%             .z, .u, .rho
%
% OUTPUT:
%   sys : struct with fields
%       .A_r, .B_r
%       .A_s, .B_s
%===========================================================

    Ns = numel(samples);

    %-------------------------------------------------------
    % Reference matrices
    %-------------------------------------------------------
    A_r = model.fun.A(z_ref, u_ref, rho_nom);
    B_r = model.fun.B_mod(z_ref, u_ref, rho_nom);

    nx = size(A_r,1);
    nu = size(B_r,2);

    %-------------------------------------------------------
    % Preallocate sample matrices
    %-------------------------------------------------------
    A_s = zeros(nx, nx, Ns);
    B_s = zeros(nx, nu, Ns);

    %-------------------------------------------------------
    % Sample matrices
    %-------------------------------------------------------
    for s = 1:Ns
        z_s   = samples(s).z;
        u_s   = samples(s).u;
        rho_s = samples(s).rho;

        A_s(:,:,s) = model.fun.A(z_s, u_s, rho_s);
        B_s(:,:,s) = model.fun.B_mod(z_s, u_s, rho_s);
    end

    %-------------------------------------------------------
    % Output struct
    %-------------------------------------------------------
    sys.A_r = A_r;
    sys.B_r = B_r;
    sys.A_s = A_s;
    sys.B_s = B_s;
end

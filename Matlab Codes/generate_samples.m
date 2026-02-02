function samples = generate_samples(z_ref, u_ref, rho_nom, uncert, N, radius, seed)
%===========================================================
% GENERATE_SAMPLES_WITH_INPUT_NUMERIC
%
% Generates N samples (z_s, rho_s, u_s) where:
%   z_s   ∈ ball of radius "radius" around z_ref
%   rho_s ∈ uncertainty set around rho_nom
%
% All outputs are numeric.
%
% INPUTS:
%   z_ref   : reference state (numeric vector)
%   u_ref   : reference input (numeric vector)
%   rho_nom : nominal parameters (numeric vector)
%   uncert  : relative uncertainty (e.g., 0.1 = ±10%)
%   N       : number of samples
%   radius  : state sampling radius
%   gamma   : LS regularization
%   seed    : RNG seed (optional)
%
% OUTPUT:
%   samples : struct array with fields
%       .z
%       .u
%       .rho
%       .e
%===========================================================

    if nargin >= 7 && ~isempty(seed)
        rng(seed);
    end

    nz = numel(z_ref);

    % Preallocate output
    samples(N,1) = struct( ...
        'z', [], ...
        'u', [], ...
        'rho', [], ...
        'e', []);

    % --- Sampling loop
    for s = 1:N
        %-------------------------------------------------------
        % 1. Sample state deviation in n-ball
        %-------------------------------------------------------
        dir = randn(nz,1);
        dir = dir / norm(dir);
        mag = radius * rand^(1/nz);
        e   = mag * dir;

        z_s = z_ref + e;

        %-------------------------------------------------------
        % 2. Sample uncertain parameters
        %-------------------------------------------------------
        rho_s = rho_nom .* (1 + uncert .* (2*rand(size(rho_nom)) - 1));

        %-------------------------------------------------------
        % 3. Store numeric sample
        %-------------------------------------------------------
        samples(s).z   = z_s;
        samples(s).u   = u_ref;
        samples(s).rho = rho_s;
        samples(s).e   = e;
    end
end
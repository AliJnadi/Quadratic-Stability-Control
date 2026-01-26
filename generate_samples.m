function samples = generate_samples(x_ref, rho_nom, uncert, N, radius, seed)
%GENERATE_SAMPLES
% Generates N samples (x_i, rho_i) with
%   x_i ∈ B_radius(x_ref)
%   rho_i ∈ admissible uncertainty set around rho_nom
%
% Sampling of state and parameters is independent.

    if nargin >= 6
        rng(seed);
    end

    n = numel(x_ref);

    % Preallocate struct array
    samples(N,1) = struct('x', [], 'e', [], 'rho', []);

    for s = 1:N
        % 1. Sample state deviation uniformly in an n-ball
        dir = randn(n,1);
        dir = dir / norm(dir);                  % random direction
        mag = radius * rand^(1/n);              % correct radial distribution
        e   = mag * dir;                        % state deviation

        x_s = x_ref + e;

        % 2. Sample uncertain parameters
        rho_s = sample_rho(rho_nom, uncert);

        % 3. Store sample
        samples(s).x   = x_s;
        samples(s).e   = e;
        samples(s).rho = rho_s;
    end
end

function rho = sample_rho(rho_nom, uncert)
% Samples rho from a hyper-rectangle centered at rho_nom

    rho = rho_nom;
    fields = fieldnames(rho_nom);

    for k = 1:numel(fields)
        f = fields{k};
        v = rho_nom.(f);
        rho.(f) = v * (1 + uncert * (2*rand - 1));
    end
end
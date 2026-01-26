function samples = generate_samples_with_input(model, z_ref, u_nom, rho_nom, uncert, N, radius, gamma, force_selector, seed)
%===========================================================
% GENERATE_SAMPLES_WITH_INPUT_NUMERIC
%
% Generates N samples (z_s, rho_s, u_s) where:
%   z_s   ∈ ball of radius "radius" around z_ref
%   rho_s ∈ uncertainty set around rho_nom
%   u_s   minimizes || f(z_s, u, rho_s) ||
%
% All outputs are numeric.
%
% INPUTS:
%   model   : symbolic model struct (fields f, B_mod, z, u_mod, params)
%   z_ref   : reference state (numeric vector)
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
%       .e
%       .rho
%       .u
%       .residual_norm
%===========================================================

    if nargin >= 10 && ~isempty(seed)
        rng(seed);
    end

    nz = numel(z_ref);
    nu_mod = length(model.u_mod);

    % Preallocate output
    samples(N,1) = struct( ...
        'z', [], ...
        'e', [], ...
        'rho', [], ...
        'u', [], ...
        'residual_norm', [] );

%     % --- Convert symbolic model to numeric functions (once)
%     f_fun = model.fun.f;
%     B_fun = model.fun.B_mod;

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
        
        options = optimoptions(...
                               'fsolve', ...
                               'Algorithm', 'levenberg-marquardt', ...
                               'Display', 'off');
                           
        [u_opt, ~, ~] = fsolve(@(u) model.fun.f(z_s, u, rho_s), u_nom, options);
        
        % Split results
        mask = logical([force_selector, force_selector]);
        u_s = u_opt(mask);
        
        u_mod = u_opt;
        u_mod(mask == 0) = 0;
        res = model.fun.f(z_s, u_mod, rho_s);

%         %-------------------------------------------------------
%         % 3. Evaluate f(z_s,0,rho_s) numerically
%         %-------------------------------------------------------
%         f0 = f_fun(z_s, zeros(length(model.u),1), rho_s);
% 
%         %-------------------------------------------------------
%         % 4. Evaluate B_mod(z_s,0,rho_s) numerically
%         %-------------------------------------------------------
%         B  = B_fun(z_s, zeros(nu_mod,1), rho_s);

%         %-------------------------------------------------------
%         % 5. Least-squares input to minimize ||f + B u||
%         %-------------------------------------------------------
%         u_s = - (B.'*B + gamma*eye(nu_mod)) \ (B.'*f0);

%         %-------------------------------------------------------
%         % 6. Residual
%         %-------------------------------------------------------
%         res = f0 + B*u_s;

        %-------------------------------------------------------
        % 7. Store numeric sample
        %-------------------------------------------------------
        samples(s).z = z_s;
        samples(s).e = e;
        samples(s).rho = rho_s;
        samples(s).u = u_s;
        samples(s).residual_norm = norm(res);
    end
end

% function samples = generate_samples_with_input(model, z_ref, rho_nom, uncert, N, radius, gamma, seed)
% %===========================================================
% % GENERATE_SAMPLES_WITH_INPUT
% %
% % Generates N samples (z_s, rho_s, u_s) where:
% %   z_s   ∈ ball of radius "radius" around z_ref
% %   rho_s ∈ uncertainty set around rho_nom
% %   u_s   minimizes || f(z_s, u, rho_s) ||
% %
% % INPUTS:
% %   model   : symbolic model struct (fields f, z, u, params, B)
% %   z_ref   : reference state
% %   rho_nom : nominal parameter struct
% %   uncert  : relative uncertainty (e.g. 0.1 = ±10%)
% %   N       : number of samples
% %   radius  : state sampling radius
% %   gamma   : LS regularization parameter
% %   seed    : RNG seed (optional)
% %
% % OUTPUT:
% %   samples : struct array with fields
% %       .z
% %       .e
% %       .rho
% %       .u
% %       .residual_norm
% %===========================================================
% 
%     if nargin >= 8 && ~isempty(seed)
%         rng(seed);
%     end
%     
%     z = model.z;
%     u = model.u;
%     u_mod = model.u_mod;
%     p = model.params;
%     
%     nz = numel(z_ref);
%     nu = length(u);
%     nu_mod = length(u_mod);
% 
%     % Preallocate
%     samples(N,1) = struct( ...
%         'z', [], ...
%         'e', [], ...
%         'rho', [], ...
%         'u', [], ...
%         'residual_norm', [] );
% 
%     for s = 1:N
%         %-------------------------------------------------------
%         % 1. Sample state deviation in an n-ball
%         %-------------------------------------------------------
%         dir = randn(nz,1);
%         dir = dir / norm(dir);
%         mag = radius * rand^(1/nz);
%         e   = mag * dir;
% 
%         z_s = z_ref + e;
% 
%         %-------------------------------------------------------
%         % 2. Sample uncertain parameters
%         %-------------------------------------------------------
%         rho_s = sample_rho(rho_nom, uncert);
% 
%     %     z_nom = double(subs(z_eq, model.params, rho_nom));
%         u_nom = zeros(nu,1);
%         u_mod = zeros(nu_mod,1);
% 
%         %-------------------------------------------------------
%         % 3. Evaluate nonlinear drift f(z_s,0,rho_s)
%         %-------------------------------------------------------
%         f0 = subs(model.f, ...
%                   [z; u; p], ...
%                   [z_s; u_nom; rho_s]);
%     %     f0 = double(subs( ...
%     %         model.f, ...
%     %         [model.z; model.u; model.params], ...
%     %         [z_s; zeros(nu,1); rho_struct_to_vector(rho_s, model.params)] ));
%         f0 = simplify(f0, 'Steps', 30);
% 
%         %-------------------------------------------------------
%         % 4. Evaluate input Jacobian B(z_s,0,rho_s)
%         %-------------------------------------------------------
%         B = subs(model.B_mod, ...
%                  [z; u_mod; p], ...
%                  [z_s; u_mod; rho_s]);
% 
%         B = simplify(B, 'Steps', 30);
%     %     B = double(subs( ...
%     %         model.B_mod, ...
%     %         [model.z; model.u; model.params], ...
%     %         [z_s; zeros(nu,1); rho_struct_to_vector(rho_s, model.params)] ));
% 
%         %-------------------------------------------------------
%         % 5. Least-squares cancellation input
%         %-------------------------------------------------------
%         u_s = - (B.'*B + gamma*eye(nu_mod)) \ (B.'*f0);
% 
%         %-------------------------------------------------------
%         % 6. Residual (diagnostic)
%         %-------------------------------------------------------
%         res = f0 + B*u_s;
% 
%         %-------------------------------------------------------
%         % 7. Store sample
%         %-------------------------------------------------------
%         samples(s).z = z_s;
%         samples(s).e = e;
%         samples(s).rho = rho_s;
%         samples(s).u = u_s;
%         samples(s).residual_norm = norm(res);
%     end
% end
% 
% function rho = sample_rho(rho_nom, uncert)
% %SAMPLE_RHO
% % Samples a parameter vector rho from a hyper-rectangle
% % centered at rho_nom with relative uncertainty 'uncert'.
%     rho = rho_nom .* (1 + uncert .* (2*rand(size(rho_nom)) - 1));
% end
% 

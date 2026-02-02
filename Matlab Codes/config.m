% config.m - system configuration

% Diagonistic flag
% 0 = silent
% 1 = high-level info
% 2 = detailed diagnostics
verbose = 1;

Q_eps = 1e-4;   % Lyapunov strictness parameter

% ODE
eps_f = 1e-12;
eps_z = 1e-8;

t_stable = pi;
Tspan = [0, 15];
opts = odeset( ...
    'RelTol',        1e-6, ...
    'AbsTol',        1e-8, ...
    'MaxStep',       0.05, ...
    'InitialStep',  1e-3, ...
    'Refine',        4 ...
);
 
N = 100;            % Number of samples used for LMI
N_t = 30;           % Number of samples used for testing
r = 0.05;           % Sampling radius
uncert = 0.1;       % Uncertainty  
seed = 1;           % Random seed
gamma = 1e-2;       % Regularization for LS

u_mag = 0.25;   % Deviation from equilibrium

row = 2; 
col = 2;

force_selector = [1 1 1 1];  % forces applied on m11

% System nominal parameters
m_val = ones(row*col,1);

l = 1;
k = 1;
mu_l = 1;

kd = 0.5;
mu_d = 0.5;

rho_nom = [ ...
    l;      % l
    k;      % k
    mu_l;   % mu_l
    kd;    % kd
    mu_d;   % mu_d
    m_val;  % masses
];
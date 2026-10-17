function [z0, report] = validate_symbolic_model(model, row, col, verbose)
%===========================================================
% VALIDATE_GRID_MODEL
%
% Algebraic (non-ODE) validation of symbolic grid model
%
% VERBOSE:
%   0 -> no output
%   1 -> pass/fail only
%   2 -> detailed diagnostic output
%===========================================================

if nargin < 4
    verbose = 1;
end

n = row;
m = col;
N = n*m;

% -------- Extract symbols --------
f = model.sym.f;
z = model.z;
u = model.u;

% -------- Construct symbolic equilibrium --------
syms l real

% Expected equilibirium points     
N = row * col;
[x0, y0] = meshgrid(l:l:col*l, l:l:row*l);
x0 = reshape(x0',N,1);
y0 = reshape(y0',N,1);

z0 = [x0; y0; zeros(2*N,1)];
u0 = zeros(size(u));
% =========================================================
% Exact equilibrium check
% =========================================================
f_eq = double(subs(f, [z; u], [z0; u0]));
eq_ok = isequal(f_eq, zeros(size(f)));

% -------- Pack report --------
report.valid = eq_ok;
report.f_at_equilibrium = f_eq;
report.z_equilibrium = z0;

% =========================================================
% Verbose output control
% =========================================================
if verbose > 0
    disp('------------------------------------------')
    disp('GRID MODEL VALIDATION')
    disp('------------------------------------------')
    disp(['Model Valid   : ', tf(eq_ok)])
end

if verbose == 2
    disp(' ')
    disp('--- Detailed diagnostics ---')

    if ~eq_ok
        disp('f(z*,u*) ≠ 0 :')
        disp(f_eq)
    else
        disp('f(z*,u*) = 0 ✓')
    end

    disp('Equilibrium state x*:')
    x0 = z0(1:N);
    disp(reshape(x0, row, col)');
    
    disp('Equilibrium state y*:')
    y0 = z0(N+1:2*N);
    disp(reshape(y0, row, col)');
end

end

% -------- Helper --------
function s = tf(b)
    if b
        s = 'TRUE';
    else
        s = 'FALSE';
    end
end
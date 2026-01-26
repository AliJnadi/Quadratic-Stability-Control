function out = f_grid_symbolic(row, col, force_selector, verbose)
%===========================================================
% Exact nonlinear symbolic model of an n×m mass–spring–damper grid
% with fixed walls aligned to rest length l.
%
% OUTPUT (struct):
%   out.f      : nonlinear symbolic vector field
%   out.z      : state vector
%   out.u      : input vector
%   out.A      : Jacobian df/dz
%   out.B      : Jacobian df/du
%   out.params : symbolic parameters
%
% Recommended usage: row = col = 2 force_selector = [1 0 0 0]
%===========================================================

if nargin < 4
    verbose = 1;
end

if nargin >= 8 && ~isempty(seed)
    rng(seed);
end

n = row;
m = col;
N = n*m;

% -------- Symbolic variables --------
syms l real
syms K mu_l Kd mu_d real
syms m_ij [N 1] real
ld = sqrt(sym(2))*l;

% states
syms x  [N 1] real
syms y  [N 1] real
syms dx [N 1] real
syms dy [N 1] real

z = [x; y; dx; dy];

% inputs (forces)
syms ux [N 1] real
syms uy [N 1] real
u = [ux; uy];

params = [l; K; mu_l; Kd; mu_d; m_ij];

% -------- Extended grid with fixed walls --------  
% walls: exact multiples of l
[xE, yE] = meshgrid(0:l:(col + 1)*l, 0:l:(row + 1)*l);
dxE = sym(zeros(n+2,m+2));
dyE = sym(zeros(n+2,m+2));

xE(2:end-1,2:end-1)  = reshape(x,n,m)';
yE(2:end-1,2:end-1)  = reshape(y,n,m)';
dxE(2:end-1,2:end-1) = reshape(dx,n,m)';
dyE(2:end-1,2:end-1) = reshape(dy,n,m)';

if verbose >= 2
    display(x);
    display(y);
    display(dx);
    display(dy);
    display(xE);
    display(yE);
    display(dxE);
    display(dyE);
end

% -------- Force balance --------
ddx = sym(zeros(N,1));
ddy = sym(zeros(N,1));
id = 1;
for i = 2:n+1
  for j = 2:m+1
    
    r1 = sqrt((xE(i, j) - xE(i-1, j-1))^2 + (yE(i, j) - yE(i-1, j-1))^2);
    r2 = sqrt((xE(i, j) - xE(i-1, j+1))^2 + (yE(i, j) - yE(i-1, j+1))^2);
    r3 = sqrt((xE(i, j) - xE(i+1, j-1))^2 + (yE(i, j) - yE(i+1, j-1))^2);
    r4 = sqrt((xE(i, j) - xE(i+1, j+1))^2 + (yE(i, j) - yE(i+1, j+1))^2);
    
    % On x direction
    % Input 
    Fx = ux(id);
    
    % Linear Spring Damper
    Fx = Fx - K*((xE(i, j)  - xE(i, j-1)) - l) - mu_l * (dxE(i, j) - dxE(i, j-1));
    Fx = Fx + K*((-xE(i, j) + xE(i, j+1)) - l) - mu_l * (dxE(i, j) - dxE(i, j+1));
    
    % Diagonal Spring Damper
    Fx = Fx - Kd*(r1 - ld)* (xE(i, j)  - xE(i - 1, j - 1))/ r1 - mu_d * (dxE(i, j) - dxE(i - 1, j - 1));
    Fx = Fx - Kd*(r2 - ld)* (xE(i, j)  - xE(i - 1, j + 1))/ r2 - mu_d * (dxE(i, j) - dxE(i - 1, j + 1));
    Fx = Fx + Kd*(r3 - ld)* (-xE(i, j) + xE(i + 1, j - 1))/ r3 - mu_d * (dxE(i, j) - dxE(i + 1, j - 1));
    Fx = Fx + Kd*(r4 - ld)* (-xE(i, j) + xE(i + 1, j + 1))/ r4 - mu_d * (dxE(i, j) - dxE(i + 1, j + 1));
            
    % On y direction
    % Input 
    Fy = uy(id);
    
    % Linear Spring Damper
    Fy = Fy - K*((yE(i, j)  - yE(i-1, j)) - l) - mu_l * (dyE(i, j) - dyE(i-1, j));
    Fy = Fy + K*((-yE(i, j) + yE(i+1, j)) - l) - mu_l * (dyE(i, j) - dyE(i+1, j));
    
    % Diagonal Spring Damper
    Fy = Fy - Kd*(r1 - ld)* (yE(i, j)  - yE(i - 1, j - 1))/ r1 - mu_d * (dyE(i, j) - dyE(i - 1, j - 1));
    Fy = Fy - Kd*(r2 - ld)* (yE(i, j)  - yE(i - 1, j + 1))/ r2 - mu_d * (dyE(i, j) - dyE(i - 1, j + 1));
    Fy = Fy + Kd*(r3 - ld)* (-yE(i, j) + yE(i + 1, j - 1))/ r3 - mu_d * (dyE(i, j) - dyE(i + 1, j - 1));
    Fy = Fy + Kd*(r4 - ld)* (-yE(i, j) + yE(i + 1, j + 1))/ r4 - mu_d * (dyE(i, j) - dyE(i + 1, j + 1));

    ddx(id) = Fx / m_ij(id);
    ddy(id) = Fy / m_ij(id);
    
    id = id + 1;
  end
end

% -------- Nonlinear state equation --------
f = [
    dx
    dy
    ddx
    ddy
];

% -------- Jacobians (exact, general) --------
A = jacobian(f, z);
B = jacobian(f, u);

% ---------- Applying Force Selector ----------
[B_mod, u_mod] = apply_force_selector(B, u, force_selector);

% Numeric function handles
f_fun = matlabFunction(f, 'Vars', {z, u, params});
A_fun = matlabFunction(A, 'Vars', {z, u, params});
B_fun = matlabFunction(B, 'Vars', {z, u, params});
B_mod_fun = matlabFunction(B_mod, 'Vars', {z, u_mod, params});

% -------- Output --------
out.sym.f      = simplify(f,'Steps',50);
out.sym.A      = simplify(A,'Steps',50);
out.sym.B      = simplify(B,'Steps',50);
out.sym.B_mod  = simplify(B_mod,'Steps',50);

out.fun.f     = f_fun;
out.fun.A     = A_fun;
out.fun.B     = B_fun;
out.fun.B_mod = B_mod_fun;

out.z      = z;
out.u      = u;
out.u_mod  = u_mod;

out.params = params;
end

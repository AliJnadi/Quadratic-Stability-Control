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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %    r1 r2 r3    
    %    r4 ** r6
    %    r7 r8 r9
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % distances
    r2 = hypot(xE(i,j)-xE(i-1,j), yE(i,j)-yE(i-1,j));
    r8 = hypot(xE(i,j)-xE(i+1,j), yE(i,j)-yE(i+1,j));
    r4 = hypot(xE(i,j)-xE(i,j-1), yE(i,j)-yE(i,j-1));
    r6 = hypot(xE(i,j)-xE(i,j+1), yE(i,j)-yE(i,j+1));

    r1 = hypot(xE(i,j)-xE(i-1,j-1), yE(i,j)-yE(i-1,j-1));
    r3 = hypot(xE(i,j)-xE(i-1,j+1), yE(i,j)-yE(i-1,j+1));
    r7 = hypot(xE(i,j)-xE(i+1,j-1), yE(i,j)-yE(i+1,j-1));
    r9 = hypot(xE(i,j)-xE(i+1,j+1), yE(i,j)-yE(i+1,j+1));
    
    % Input 
    Fx = ux(id);
    Fy = uy(id);
    
     % axial springs
    Fx = Fx ...
        - K*(r4-l)*(xE(i,j)-xE(i,j-1))/(r4 + 1e-6) - mu_l*(dxE(i,j)-dxE(i,j-1)) ...
        - K*(r6-l)*(xE(i,j)-xE(i,j+1))/(r6 + 1e-6)  - mu_l*(dxE(i,j)-dxE(i,j+1));

    Fy = Fy ...
        - K*(r2-l)*(yE(i,j)-yE(i-1,j))/(r2 + 1e-6) - mu_l*(dyE(i,j)-dyE(i-1,j)) ...
        - K*(r8-l)*(yE(i,j)-yE(i+1,j))/(r8 + 1e-6) - mu_l*(dyE(i,j)-dyE(i+1,j));

    % diagonal springs
    Fx = Fx ...
        - Kd*(r1-ld)*(xE(i,j)-xE(i-1,j-1))/(r1 + 1e-6) - mu_d*(dxE(i,j)-dxE(i-1,j-1)) ...
        - Kd*(r3-ld)*(xE(i,j)-xE(i-1,j+1))/(r3 + 1e-6) - mu_d*(dxE(i,j)-dxE(i-1,j+1)) ...
        - Kd*(r7-ld)*(xE(i,j)-xE(i+1,j-1))/(r7 + 1e-6) - mu_d*(dxE(i,j)-dxE(i+1,j-1)) ...
        - Kd*(r9-ld)*(xE(i,j)-xE(i+1,j+1))/(r9 + 1e-6) - mu_d*(dxE(i,j)-dxE(i+1,j+1));

    Fy = Fy ...
        - Kd*(r1-ld)*(yE(i,j)-yE(i-1,j-1))/(r1 + 1e-6) - mu_d*(dyE(i,j)-dyE(i-1,j-1)) ...
        - Kd*(r3-ld)*(yE(i,j)-yE(i-1,j+1))/(r3 + 1e-6) - mu_d*(dyE(i,j)-dyE(i-1,j+1)) ...
        - Kd*(r7-ld)*(yE(i,j)-yE(i+1,j-1))/(r7 + 1e-6) - mu_d*(dyE(i,j)-dyE(i+1,j-1)) ...
        - Kd*(r9-ld)*(yE(i,j)-yE(i+1,j+1))/(r9 + 1e-6) - mu_d*(dyE(i,j)-dyE(i+1,j+1));

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
out.sym.f      = f;
out.sym.A      = A;
out.sym.B      = B;
out.sym.B_mod  = B_mod;

out.fun.f     = f_fun;
out.fun.A     = A_fun;
out.fun.B     = B_fun;
out.fun.B_mod = B_mod_fun;

out.z      = z;
out.u      = u;
out.u_mod  = u_mod;

out.params = params;
end

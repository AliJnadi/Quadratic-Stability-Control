function out = f_grid_symbolic_v1(row, col, force_selector)
%===========================================================
% Exact nonlinear symbolic model of an n×m mass–spring–damper grid
% with fixed walls aligned to rest length l.
% Force selector determines which masses are actuated (1 = active).
%
% OUTPUT (struct):
%   out.f        : nonlinear symbolic vector field
%   out.z        : state vector
%   out.u        : input vector (full)
%   out.A        : Jacobian df/dz
%   out.B        : Jacobian df/du (full)
%   out.B_mod    : Jacobian df/du (only active forces)
%   out.params   : symbolic parameters
%
%===========================================================

n = row; m = col; N = n*m;

%% Symbolic variables
syms l real
syms K mu_l Kd mu_d real
syms m_ij [N 1] real

% states
syms x [N 1] real
syms y [N 1] real
syms dx [N 1] real
syms dy [N 1] real
z = [x; y; dx; dy];

% inputs (forces)
syms ux [N 1] real
syms uy [N 1] real
u = [ux; uy];

params = [l; K; mu_l; Kd; mu_d; m_ij];

%% Indexing helper
idx = @(i,j) (j-1)*n + i;

%% Extended grid with walls (symbolic)
xE  = sym(zeros(n+2,m+2));
yE  = sym(zeros(n+2,m+2));
dxE = sym(zeros(n+2,m+2));
dyE = sym(zeros(n+2,m+2));

xE(2:end-1,2:end-1)  = reshape(x,n,m);
yE(2:end-1,2:end-1)  = reshape(y,n,m);
dxE(2:end-1,2:end-1) = reshape(dx,n,m);
dyE(2:end-1,2:end-1) = reshape(dy,n,m);

% Wall positions: exact multiples of l
for j = 1:m+2
    xE(:,j) = (j-1)*l;
end
for i = 1:n+2
    yE(i,:) = (n+2-i)*l;
end

%% Neighbor list: [di,dj,K,mu,l0]
ld = sqrt(2)*l;
nbr = [
 0  1  K  mu_l  l
 0 -1  K  mu_l  l
 1  0  K  mu_l  l
-1  0  K  mu_l  l
 1  1  Kd mu_d ld
 1 -1  Kd mu_d ld
-1  1  Kd mu_d ld
-1 -1  Kd mu_d ld
];

%% Spring-damper function
% spring = @(xi,yi,dxi,dyi,xk,yk,dxk,dyk,Ks,cs,ls) ...
%     Ks*((sqrt((xk-xi)^2 + (yk-yi)^2) - ls)/sqrt((xk-xi)^2 + (yk-yi)^2))*[xk-xi; yk-yi] + ...
%     cs*(((dxk-dxi)*(xk-xi) + (dyk-dyi)*(yk-yi))/((xk-xi)^2 + (yk-yi)^2))*[xk-xi; yk-yi];
spring = @(xi,yi,dxi,dyi,xk,yk,dxk,dyk,Ks,cs,ls) ...
    Ks * ((sqrt((xk-xi)^2 + (yk-yi)^2) - ls)/sqrt((xk-xi)^2 + (yk-yi)^2)) * [xk-xi; yk-yi] + ...
    cs * (((dxk-dxi)*(xk-xi) + (dyk-dyi)*(yk-yi)) / ((xk-xi)^2 + (yk-yi)^2)) * [xk-xi; yk-yi];


%% Initialize acceleration
ddx = sym(zeros(N,1));
ddy = sym(zeros(N,1));

%% Force accumulation
for j = 1:m
    for i = 1:n
        id = idx(i,j);
        Fx = ux(id);
        Fy = uy(id);
        
        for k = 1:8
            ii = i + nbr(k,1);
            jj = j + nbr(k,2);
            
            Ks = nbr(k,3);
            cs = nbr(k,4);
            ls = nbr(k,5);
            
            % skip neighbors outside grid (walls will be handled)
            if ii>=1 && ii<=n && jj>=1 && jj<=m
                id2 = idx(ii,jj);
                
                F = spring( ...
                    xE(i+1,j+1), yE(i+1,j+1), dxE(i+1,j+1), dyE(i+1,j+1), ...
                    xE(ii+1,jj+1), yE(ii+1,jj+1), dxE(ii+1,jj+1), dyE(ii+1,jj+1), ...
                    Ks, cs, ls);
                
                % Newton's 3rd law
                ddx(id)  = ddx(id)  + F(1)/m_ij(id);
                ddy(id)  = ddy(id)  + F(2)/m_ij(id);
                ddx(id2) = ddx(id2) - F(1)/m_ij(id2);
                ddy(id2) = ddy(id2) - F(2)/m_ij(id2);
            else
                % Mass-wall spring (wall is fixed)
                F = spring( ...
                    xE(i+1,j+1), yE(i+1,j+1), dxE(i+1,j+1), dyE(i+1,j+1), ...
                    xE(ii+1,jj+1), yE(ii+1,jj+1), 0, 0, ...
                    Ks, cs, ls);
                
                ddx(id) = ddx(id) + F(1)/m_ij(id);
                ddy(id) = ddy(id) + F(2)/m_ij(id);
            end
        end
    end
end

%% Nonlinear state equation
f = [dx; dy; ddx; ddy];

%% Jacobians
A = jacobian(f, z);
B = jacobian(f, u);

%% Apply force selector (2 DOF per selected mass)
B_mod = apply_force_selector(B, force_selector);

%% Output
out.f      = simplify(f,'Steps',50);
out.z      = z;
out.u      = u;
out.A      = simplify(A,'Steps',50);
out.B      = simplify(B,'Steps',50);
out.B_mod  = simplify(B_mod,'Steps',50);
out.params = params;

end
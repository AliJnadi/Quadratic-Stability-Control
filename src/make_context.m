function ctx=make_context(fs,grid)
if nargin<1,fs=ones(1,4);end
if nargin<2,grid=2;end
config;row=grid;col=grid;rho_nom=[l;k;mu_l;kd;mu_d;ones(grid^2,1)];
model=f_grid_symbolic(row,col,fs,0);
[eq,rep]=validate_symbolic_model(model,row,col,0);assert(rep.valid);
[zref,uref,ufull,refinfo]=generate_operating_point(model,eq,rho_nom,u_mag,opts,eps_f,eps_z,fs,seed,0,row,col,l);
assert(refinfo.residual_norm<1e-9);
ctx=struct('model',model,'fs',fs,'zref',zref,'uref',uref,'ufull',ufull,'rho_nom',rho_nom,'reference',refinfo,'grid',grid);
end

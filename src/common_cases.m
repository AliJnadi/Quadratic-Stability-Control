function cases=common_cases(ctx,cfg)
u=generate_samples(ctx.zref,ctx.uref,ctx.rho_nom,cfg.uncertainty,cfg.uniform_count,cfg.radius,cfg.validation_seed,1);
s=generate_samples(ctx.zref,ctx.uref,ctx.rho_nom,cfg.uncertainty,cfg.stress_count,cfg.radius,cfg.validation_seed+1,0);
all=[u;s];cases=struct('id',{},'set',{},'scenario',{},'z0',{},'rho',{},'disturbed',{});
for scenario=1:3
 for j=1:numel(all)
  label='uniform';if j>numel(u),label='stress';end
  rho=all(j).rho;z0=all(j).z;disturbed=false;
  if scenario==1,rho=ctx.rho_nom;name='nominal_regulation';
  elseif scenario==2,name='uncertain_regulation';
  else,name='uncertain_disturbance';z0=ctx.zref;disturbed=true;end
  cases(end+1)=struct('id',numel(cases)+1,'set',label,'scenario',name,'z0',z0,'rho',rho,'disturbed',disturbed);
 end
end
cases(end+1)=struct('id',numel(cases)+1,'set','unique_nominal','scenario','nominal_disturbance','z0',ctx.zref,'rho',ctx.rho_nom,'disturbed',true);
end

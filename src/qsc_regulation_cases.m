function cases = qsc_regulation_cases(ctx,samples)
cases = struct('id',{},'set',{},'scenario',{},'z0',{},'rho',{},'disturbed',{});
for nominal = [true false]
 for j = 1:numel(samples)
  rho = samples(j).rho; scenario = 'uncertain_regulation';
  if nominal, rho = ctx.rho_nom; scenario = 'nominal_regulation'; end
  cases(end+1) = struct('id',numel(cases)+1,'set','boundary_corner', ...
   'scenario',scenario,'z0',samples(j).z,'rho',rho,'disturbed',false); %#ok<AGROW>
 end
end
end

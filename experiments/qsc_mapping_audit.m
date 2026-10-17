function qsc_mapping_audit(cfg,folder)
rows = {};
for mask = cfg.mapping_masks
 ctx = qsc_context(mask); model = ctx.model; fs = ctx.fs;
 [S,indices] = actuator_mapping(fs);
 assert(isequal(indices,[find(fs(:));4+find(fs(:))]));
 rng(9900+mask);
 for trial = 1:4
  z = ctx.zref+0.01*randn(16,1);
  rho = ctx.rho_nom.*(1+0.1*(2*rand(9,1)-1));
  u = 0.1*randn(size(S,2),1); h = 1e-6;
  B = model.fun.B_mod(z,u,rho); A = model.fun.A(z,u,rho);
  J = zeros(size(B)); JA = zeros(16);
  for j = 1:numel(u)
   du = zeros(size(u)); du(j) = h;
   J(:,j) = (model.fun.f(z,S*(u+du),rho)-model.fun.f(z,S*(u-du),rho))/(2*h);
  end
  for j = 1:16
   dz = zeros(16,1); dz(j) = h;
   JA(:,j) = (model.fun.f(z+dz,S*u,rho)-model.fun.f(z-dz,S*u,rho))/(2*h);
  end
  be = max(abs(J(:)-B(:))); ae = max(abs(JA(:)-A(:)));
  assert(norm(B-model.fun.B(z,S*u,rho)*S,'fro')<1e-12 && be<1e-8 && ae<1e-7);
  rows(end+1,:) = {mask,numel(u),trial,be,ae}; %#ok<AGROW>
 end
 % Independent sign check against MATLAB's Riccati gain.
 A = model.fun.A(ctx.zref,ctx.uref,ctx.rho_nom);
 B = model.fun.B_mod(ctx.zref,ctx.uref,ctx.rho_nom);
 [native,P] = lqr(A,B,eye(16),eye(size(B,2))); K = -native;
 assert(max(real(eig(A+B*K)))<0);
 residual = A'*P+P*A-P*B*B'*P+eye(16);
 assert(norm(residual,'fro')<1e-7);
 % Check A+B*K against the nonlinear closed-loop derivative.
 closed = @(x) model.fun.f(x,expand_input(ctx.uref+K*(x-ctx.zref),fs),ctx.rho_nom);
 J = zeros(16);
 for j = 1:16
  dz = zeros(16,1); dz(j) = h;
  J(:,j) = (closed(ctx.zref+dz)-closed(ctx.zref-dz))/(2*h);
 end
 assert(max(max(abs(J-A-B*K)))<1e-7);
 fprintf('Verified mapping and LQR sign, mask %d\n',mask);
end
[S,indices] = actuator_mapping(zeros(1,4));
assert(isempty(indices) && isequal(size(S),[8 0]));
assert(isequal(expand_input(zeros(0,1),zeros(1,4)),zeros(8,1)));
T = cell2table(rows,'VariableNames',{'mask','inputs','trial','B_error','A_error'});
writetable(T,fullfile(folder,'mapping_audit.csv'));
end

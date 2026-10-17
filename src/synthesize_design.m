function design=synthesize_design(ctx,samples,cfg,kind)
sys=compute_system_matrices(ctx.model,ctx.zref,ctx.uref,ctx.rho_nom,samples);
n=size(sys.A_r,1);m=size(sys.B_r,2);p=n/2;
C=[eye(p)/cfg.position_scale zeros(p);zeros(m,n)];D=[zeros(p,m);eye(m)/cfg.force_scale];
E=zeros(n,p,numel(samples)+1);E(:,:,1)=ctx.model.fun.B(ctx.zref,ctx.ufull,ctx.rho_nom)*cfg.disturbance_scale;
for s=1:numel(samples),E(:,:,s+1)=ctx.model.fun.B(samples(s).z,ctx.ufull,samples(s).rho)*cfg.disturbance_scale;end
design=struct('kind',kind,'status','NotRun','accepted',false,'elapsed',NaN,'K',[],'X',[],'Y',[],'eta',NaN,'gamma',NaN,'error','');
% Preserve the complete attempted problem even if the solve fails early.
design.sys=sys;design.C=C;design.D=D;design.E=E;design.samples=samples;design.cfg=cfg;
if any(strcmp(kind,{'QSC','Hinf'})),design.solver=cfg.active_solver;else,design.solver='MATLAB lqr (nominal Riccati)';end
tstart=tic;
try
 if strcmp(kind,'QSC')
  [K,X,info]=solve_quadratic_stability(sys,cfg.lmi_margin*eye(n),0,'QSC');
  status=info.cvx_status; elapsed=info.elapsed_time;
  if ~info.feasible,design.status=status;design.elapsed=elapsed;return;end
  Y=K*X;eta=NaN;gain_native=K;
 elseif strcmp(kind,'Hinf')
  A=cat(3,sys.A_r,sys.A_s);B=cat(3,sys.B_r,sys.B_s);d=p;q=size(C,1);cvx_clear;
  cvx_begin sdp quiet
   variable X(n,n) symmetric
   variable Y(m,n)
   variable eta nonnegative
   minimize(eta)
   X>=cfg.X_floor*eye(n);
   [cfg.Y_bound*eye(m) Y;Y' cfg.Y_bound*eye(n)]>=0;
   for s=1:size(A,3)
    H=A(:,:,s)*X+B(:,:,s)*Y;
    [H+H' E(:,:,s) (C*X+D*Y)';E(:,:,s)' -eta*eye(d) zeros(d,q);C*X+D*Y zeros(q,d) -eye(q)]<=-cfg.lmi_margin*eye(n+d+q);
   end
  cvx_end
  status=cvx_status;elapsed=toc(tstart);
  if ~any(strcmp(status,{'Solved','Inaccurate/Solved'})),design.status=status;design.elapsed=elapsed;return;end
  K=Y/X;gain_native=K;
 else
  if strcmp(kind,'LQR'),Q=eye(n);R=eye(m);else,Q=C'*C;R=D'*D;end
  % C'D=0 here, so there is no state/input cross term.
  [gain_native,P,poles]=lqr(sys.A_r,sys.B_r,Q,R);
  K=-gain_native;X=P\eye(n);Y=K*X;eta=NaN;status='Solved';elapsed=toc(tstart);
  design.Q=Q;design.R=R;design.nominal_poles=poles;
 end
 cert=certify_design(sys,K,X,Y,eta,C,D,E,cfg,kind);
 design.K=K;design.native_gain=gain_native;design.X=X;design.Y=Y;design.eta=eta;design.gamma=sqrt(eta);
 design.status=status;design.elapsed=elapsed;design.cert=cert;design.accepted=cert.accepted;
 design.sys=sys;design.C=C;design.D=D;design.E=E;design.samples=samples;design.cfg=cfg;
catch ex
 design.status='Exception';design.error=getReport(ex,'extended','hyperlinks','off');design.elapsed=toc(tstart);cvx_clear;
end
end

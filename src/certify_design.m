function cert=certify_design(sys,K,X,Y,eta,C,D,E,cfg,kind)
A=cat(3,sys.A_r,sys.A_s); B=cat(3,sys.B_r,sys.B_s);
n=size(A,1);m=size(B,2);ns=size(A,3);
cert.X_eigenvalues=eig((X+X')/2);cert.min_X=min(cert.X_eigenvalues);cert.max_X=max(cert.X_eigenvalues);
cert.poles=complex(zeros(n,ns));cert.stability_residual=zeros(ns,1);cert.hinf_residual=nan(ns,1);
cert.Y_norm=norm(Y,2);cert.regularizer_min_eig=min(eig([cfg.Y_bound*eye(m) Y;Y' cfg.Y_bound*eye(n)]));
for s=1:ns
 cert.poles(:,s)=eig(A(:,:,s)+B(:,:,s)*K);
 H=A(:,:,s)*X+B(:,:,s)*Y;
 cert.stability_residual(s)=max(eig((H+H')+cfg.lmi_margin*eye(n)));
 if strcmp(kind,'Hinf')
  d=size(E,2);q=size(C,1);
  M=[H+H' E(:,:,s) (C*X+D*Y)';E(:,:,s)' -eta*eye(d) zeros(d,q);C*X+D*Y zeros(q,d) -eye(q)];
  cert.hinf_residual(s)=max(eig((M+M')/2+cfg.lmi_margin*eye(size(M))));
 end
end
cert.max_pole=max(real(cert.poles(:)));
cert.max_stability_residual=max(cert.stability_residual);
cert.max_hinf_residual=max(cert.hinf_residual);
cert.accepted=all(isfinite(K(:)))&&cert.min_X>0&&cert.max_pole<-cfg.pole_tolerance;
if any(strcmp(kind,{'QSC','Hinf'}))
 cert.accepted=cert.accepted&&cert.min_X>=cfg.X_floor-cfg.X_tolerance&&cert.regularizer_min_eig>=-cfg.residual_tolerance;
 if strcmp(kind,'Hinf'), residual=cert.max_hinf_residual;else,residual=cert.max_stability_residual;end
 cert.accepted=cert.accepted&&residual<=cfg.residual_tolerance;
end
cert.imposed_margin=cfg.lmi_margin;cert.numerical_residual_tolerance=cfg.residual_tolerance;
end

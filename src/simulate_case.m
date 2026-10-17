function result=simulate_case(ctx,design,test,cfg,limit,horizon,tight,legacy)
if nargin<5,limit=Inf;end
if nargin<6,horizon=cfg.horizon;end
if nargin<7,tight=false;end
if nargin<8,legacy=false;end
p=numel(ctx.zref)/2;K=design.K;
rt=cfg.reltol;at=cfg.abstol;if tight,rt=cfg.tight_reltol;at=cfg.tight_abstol;end
opts=odeset('RelTol',rt,'AbsTol',at,'MaxStep',cfg.maxstep,'Events',@escape_event);
times=(0:cfg.dt:horizon)';
if legacy
 times=[0 horizon];opts=odeset('RelTol',1e-8,'AbsTol',1e-10,'Events',@legacy_event);
end
b=ctx.model.fun.f(ctx.zref,ctx.ufull,test.rho);
m=struct('case_id',test.id,'set',test.set,'scenario',test.scenario,'controller',design.kind,'limit_N',limit,'horizon_s',horizon,'tight',tight,'complete',false,'failure','',...
 'position_MSE_m2',NaN,'final_position_error_m',NaN,'final_velocity_error_mps',NaN,'legacy_mixed_state_MSE',NaN,'legacy_mixed_final_norm',NaN,...
 'incremental_force_RMS_N',NaN,'commanded_force_RMS_N',NaN,'applied_force_RMS_N',NaN,'peak_commanded_component_N',NaN,'peak_applied_component_N',NaN,...
 'saturation_duration_s',NaN,'saturation_fraction',NaN,'weighted_L2_ratio',NaN,'weighted_output_energy_s',NaN,'disturbance_energy_N2s',NaN,...
 'reference_velocity_norm_mps',norm(b(1:p)),'reference_acceleration_norm_mps2',norm(b(p+1:end)),...
 'terminal_velocity_norm_mps',NaN,'terminal_acceleration_norm_mps2',NaN,'tail_velocity_norm_mps',NaN,'tail_position_range_m',NaN,'settled',false,'t_end_s',NaN);
result=struct('metrics',m,'b_rho',b,'rho',test.rho,'zref',ctx.zref,'uref',ctx.uref,'z0',test.z0,'trace',[]);
try
 [t,z,te,~,ie]=ode45(@rhs,times,test.z0,opts);
 completed=abs(t(end)-horizon)<1e-8&&all(isfinite(z(:)));
 if legacy&&~isempty(te)&&ie(end)==1,completed=all(isfinite(z(:)));end
 m.complete=completed;m.t_end_s=t(end);
 if ~completed,m.failure='Integration incomplete, nonfinite, or escape threshold';end
 e=z-ctx.zref';v=e*K';command=v+ctx.uref';applied=max(-limit,min(limit,command));
 w=zeros(numel(t),p);
 if test.disturbed
  active=t>=2&t<=12;phase=t(active)-2;envelope=sin(pi*phase/10).^2;
  w(active,1)=0.1*envelope.*sin(2*pi*0.5*phase);w(active,p/2+1)=0.1*envelope.*cos(2*pi*0.5*phase);
 end
 % Performance effort uses the actual incremental applied force under clipping.
 vp=applied-ctx.uref';zp=e*design.C'+vp*design.D';
 outenergy=trapz(t,sum(zp.^2,2));physical_wenergy=trapz(t,sum(w.^2,2));wenergy=physical_wenergy/cfg.disturbance_scale^2;
 m.position_MSE_m2=trapz(t,sum(e(:,1:p).^2,2))/t(end);
 m.final_position_error_m=norm(e(end,1:p));m.final_velocity_error_mps=norm(e(end,p+1:end));
 m.legacy_mixed_state_MSE=mean(sum(e.^2,2));m.legacy_mixed_final_norm=norm(e(end,:));
 m.incremental_force_RMS_N=sqrt(trapz(t,sum(v.^2,2))/t(end));
 m.commanded_force_RMS_N=sqrt(trapz(t,sum(command.^2,2))/t(end));m.applied_force_RMS_N=sqrt(trapz(t,sum(applied.^2,2))/t(end));
 m.peak_commanded_component_N=max(abs(command(:)));m.peak_applied_component_N=max(abs(applied(:)));
 saturated=abs(command)>limit;duration_components=trapz(t,double(saturated));
 m.saturation_duration_s=trapz(t,double(any(saturated,2)));m.saturation_fraction=m.saturation_duration_s/t(end);
 m.weighted_output_energy_s=outenergy;m.disturbance_energy_N2s=physical_wenergy;
 if wenergy>0,m.weighted_L2_ratio=sqrt(outenergy/wenergy);end
 final_derivative=rhs(t(end),z(end,:)');tail=t>=0.9*t(end);
 m.terminal_velocity_norm_mps=norm(final_derivative(1:p));m.terminal_acceleration_norm_mps2=norm(final_derivative(p+1:end));
 m.tail_velocity_norm_mps=max(sqrt(sum(z(tail,p+1:end).^2,2)));
 m.tail_position_range_m=max(max(z(tail,1:p),[],1)-min(z(tail,1:p),[],1));
 m.settled=completed&&m.terminal_velocity_norm_mps<=cfg.settle_velocity&&m.terminal_acceleration_norm_mps2<=cfg.settle_acceleration&&m.tail_velocity_norm_mps<=cfg.settle_velocity&&m.tail_position_range_m<=cfg.settle_position_range;
 result.trace=struct('t',t,'z',z,'commanded',command,'applied',applied,'external_disturbance',w,'saturated',saturated,'component_saturation_duration_s',duration_components,'terminal_derivative',final_derivative);
 if isfield(cfg,'trace_stride')&&cfg.trace_stride>1
  indices=unique([1:cfg.trace_stride:numel(t) numel(t)]);
  fields={'t','z','commanded','applied','external_disturbance','saturated'};
  for field=fields,result.trace.(field{1})=result.trace.(field{1})(indices,:);end
  result.trace.quadrature_dt_s=cfg.dt;result.trace.saved_trace_dt_s=cfg.dt*cfg.trace_stride;
 end
catch ex
 m.failure=getReport(ex,'basic','hyperlinks','off');
end
if ~isfield(cfg,'save_traces') || ~cfg.save_traces
 result.trace=[];
end
result.metrics=m;
 function dz=rhs(t,z)
  force_command=ctx.uref+K*(z-ctx.zref);force_applied=max(-limit,min(limit,force_command));
  dz=ctx.model.fun.f(z,ctx.model.input_map*force_applied+input_disturbance(t,p,test.disturbed),test.rho);
 end
 function [value,terminal,direction]=escape_event(~,z)
  value=min(100-norm(z(1:p)-ctx.zref(1:p)),1000-norm(z(p+1:end)));
  terminal=1;direction=-1;
 end
 function [value,terminal,direction]=legacy_event(t,z)
  [escape,~,~]=escape_event(t,z);
  value=[norm(rhs(t,z))-1e-8;escape];terminal=[1;1];direction=[-1;-1];
 end
end

function result = qsc_legacy_trial(ctx,design,test,cfg,horizon)
% Exact adaptive-time/event law of the corrected original sweeps.
% Physical position metrics use trapezoidal time integration. Legacy metrics
% retain mixed position/velocity units and are labelled separately.
p = numel(ctx.zref)/2;
result = struct('case_id',test.id,'set',test.set,'scenario',test.scenario, ...
 'controller',design.kind,'seed',design.seed,'complete',false,'failure','', ...
 'position_MSE_m2',NaN,'final_position_error_m',NaN, ...
 'legacy_mixed_state_MSE',NaN,'legacy_mixed_final_norm',NaN, ...
 'incremental_force_RMS_N',NaN,'peak_commanded_component_N',NaN, ...
 't_end_s',NaN,'event_converged',false,'terminal_derivative_norm',NaN, ...
 'reference_acceleration_norm_mps2',NaN);
try
 options = odeset('RelTol',1e-8,'AbsTol',1e-10,'Events',@convergence);
 [t,z,te] = ode45(@rhs,[0 horizon],test.z0,options);
 e = z-ctx.zref'; v = e*design.K'; command = v+ctx.uref';
 result.complete = all(isfinite(z(:))) && (abs(t(end)-horizon)<1e-8 || ~isempty(te));
 result.position_MSE_m2 = trapz(t,sum(e(:,1:p).^2,2))/t(end);
 result.final_position_error_m = norm(e(end,1:p));
 result.legacy_mixed_state_MSE = mean(sum(e.^2,2));
 result.legacy_mixed_final_norm = norm(e(end,:));
 result.incremental_force_RMS_N = sqrt(trapz(t,sum(v.^2,2))/t(end));
 result.peak_commanded_component_N = max(abs(command(:)));
 result.t_end_s = t(end); result.event_converged = ~isempty(te);
 result.terminal_derivative_norm = norm(rhs(t(end),z(end,:)'));
 drift = ctx.model.fun.f(ctx.zref,ctx.ufull,test.rho);
 result.reference_acceleration_norm_mps2 = norm(drift(p+1:end));
 if ~result.complete, result.failure = 'Incomplete or nonfinite integration'; end
catch exception
 result.failure = getReport(exception,'basic','hyperlinks','off');
end
 function dz = rhs(~,z)
  dz = ctx.model.fun.f(z,expand_input(ctx.uref+design.K*(z-ctx.zref),ctx.fs),test.rho);
 end
 function [value,terminal,direction] = convergence(t,z)
  value = norm(rhs(t,z))-cfg.event_threshold; terminal = 1; direction = -1;
 end
end

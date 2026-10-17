function qsc_run_stage(stage,cfg,folder)
% Stage adapters preserve published protocols; all paths are project-relative.
if strcmp(stage,'mapping'), qsc_mapping_audit(cfg,folder); return; end
switch stage
 case 'sweeps', run_sweeps(cfg,folder,false);
 case 'diagnostic', run_sweeps(cfg,folder,true);
 case 'comparison', run_comparison(cfg,folder);
 case 'sampling', run_sampling(cfg,folder);
 case 'saturation', run_saturation(cfg,folder);
 case {'settling','force_resolution','saturation_settling'}
  run_repeats(stage,cfg,folder);
 case 'radius_settling', run_radius_repeats(cfg,folder);
 case 'stabilizability', run_stabilizability(cfg,folder);
 case 'bounded', run_bounded(cfg,folder);
 case {'margin','scaling'}, run_margins(stage,cfg,folder);
 otherwise, error('QSC:stage','Unknown stage %s',stage);
end
end

function run_sweeps(cfg,folder,diagnostic)
jobs = struct('key',{},'sweep',{},'value',{},'mask',{},'N',{},'radius',{},'U',{});
if diagnostic
 for U = cfg.diagnostic_uncertainties
  jobs(end+1) = job(sprintf('U%03d',round(1000*U)),'diagnostic',U,15,cfg.N,cfg.radius,U); %#ok<AGROW>
 end
else
 for N = cfg.N_values
  jobs(end+1) = job(sprintf('N%d',N),'N',N,15,N,cfg.radius,cfg.uncertainty); %#ok<AGROW>
 end
 for r = cfg.radius_values
  jobs(end+1) = job(sprintf('R%03d',round(1000*r)),'R',r,15,cfg.N,r,cfg.uncertainty); %#ok<AGROW>
 end
 for U = cfg.uncertainty_values
  jobs(end+1) = job(sprintf('U%03d',round(1000*U)),'U',U,15,cfg.N,cfg.radius,U); %#ok<AGROW>
 end
 for mask = cfg.masks
  jobs(end+1) = job(sprintf('F%02d',mask),'F',sum(dec2bin(mask,4)=='1'),mask,cfg.N,cfg.radius,cfg.uncertainty); %#ok<AGROW>
 end
end
rows = {}; audits = {};
for j = 1:numel(jobs)
 J = jobs(j); local = cfg; local.N = J.N; local.radius = J.radius;
 local.uncertainty = J.U; local.active_solver = cfg.original_solver;
 ctx = qsc_context(J.mask); sub = fullfile(folder,J.key);
 if ~exist(sub,'dir'), mkdir(sub); end
 validation = generate_samples(ctx.zref,ctx.uref,ctx.rho_nom,J.U, ...
  cfg.original_validation_count,J.radius,cfg.original_validation_seed,0);
 cases = qsc_regulation_cases(ctx,validation);
 qsc_atomic_save(fullfile(sub,'validation.mat'),struct('samples',validation,'cases',cases));
 seeds = cfg.original_design_seed; if diagnostic, seeds = cfg.seeds; end
 for kind = {'QSC','LQR'}
  design_seeds = seeds; if strcmp(kind{1},'LQR'), design_seeds = 0; end
  for seed = design_seeds
   samples = generate_samples(ctx.zref,ctx.uref,ctx.rho_nom,J.U,J.N,J.radius,max(seed,2),1);
   design = qsc_get_design(fullfile(sub,sprintf('%s_seed%d_design.mat',kind{1},seed)),ctx,samples,local,kind{1},seed);
   audits(end+1,:) = {J.key,kind{1},seed,design.status,design.accepted,design.elapsed}; %#ok<AGROW>
   if ~design.accepted, continue; end
   T = qsc_trials(ctx,design,cases,local,sub,sprintf('%s_seed%d',kind{1},seed),Inf,cfg.original_horizon,false,true);
   for scenario = {'nominal_regulation','uncertain_regulation'}
    v = T(strcmp(T.scenario,scenario{1}),:); ok = v.complete;
    rows(end+1,:) = {J.sweep,J.value,J.mask,J.N,J.radius,J.U,kind{1},seed,scenario{1}, ...
     height(v),sum(ok),mean(v.position_MSE_m2(ok)),std(v.position_MSE_m2(ok)), ...
     mean(v.final_position_error_m(ok)),std(v.final_position_error_m(ok)), ...
     mean(v.incremental_force_RMS_N(ok)),mean(v.legacy_mixed_state_MSE(ok)), ...
     mean(v.legacy_mixed_final_norm(ok)),design.elapsed}; %#ok<AGROW>
   end
   write_sweep(rows,folder); write_design_audit(audits,folder);
  end
 end
end
write_sweep(rows,folder); write_design_audit(audits,folder);
end
function j = job(key,sweep,value,mask,N,r,U)
j = struct('key',key,'sweep',sweep,'value',value,'mask',mask,'N',N,'radius',r,'U',U);
end
function write_sweep(rows,folder)
if isempty(rows), return; end
T = cell2table(rows,'VariableNames',{'sweep','value','mask','N','radius','uncertainty', ...
 'controller','seed','scenario','trials','completed','position_MSE_mean','position_MSE_std', ...
 'final_position_error_mean','final_position_error_std','incremental_force_RMS_mean', ...
 'legacy_mixed_MSE_mean','legacy_mixed_final_mean','synthesis_time_s'});
writetable(T,fullfile(folder,'summary.csv'));
end
function write_design_audit(rows,folder)
if isempty(rows), return; end
T = cell2table(rows,'VariableNames',{'key','controller','seed','status','accepted','synthesis_time_s'});
writetable(T,fullfile(folder,'controller_audit.csv'));
end

function run_comparison(cfg,folder)
ctx = qsc_context; cases = common_cases(ctx,cfg);
qsc_atomic_save(fullfile(folder,'cases.mat'),struct('cases',cases,'reference',rmfield(ctx,'model')));
qsc_case_provenance(cases,folder);
local = cfg; local.active_solver = cfg.comparison_solver; audits = {};
for kind = {'QSC','Hinf','LQR','LQR_weighted'}
 seeds = cfg.seeds; if strncmp(kind{1},'LQR',3), seeds = 0; end
 for seed = seeds
  samples = generate_samples(ctx.zref,ctx.uref,ctx.rho_nom,cfg.uncertainty,cfg.N,cfg.radius,max(seed,2),1);
  key = sprintf('%s_seed%d',kind{1},seed);
  design = qsc_get_design(fullfile(folder,[key '_design.mat']),ctx,samples,local,kind{1},seed);
  audits(end+1,:) = {key,kind{1},seed,design.status,design.accepted,design.elapsed}; %#ok<AGROW>
  if design.accepted
   qsc_trials(ctx,design,cases,local,folder,key);
   heldout(ctx,design,cases,local,folder,key);
  end
  write_design_audit(audits,folder);
 end
end
end
function [ctx,cases,files] = comparison_inputs(cfg)
ctx = qsc_context; folder = fullfile(cfg.output_dir,'comparison');
assert(exist(fullfile(folder,'completed.mat'),'file') == 2,'Run comparison stage first');
s = load(fullfile(folder,'cases.mat'),'cases'); cases = s.cases;
files = dir(fullfile(folder,'*_design.mat'));
end
function heldout(ctx,design,cases,cfg,folder,key)
cases = cases(strcmp({cases.scenario},'uncertain_regulation'));
samples = struct('z',{},'u',{},'rho',{},'e',{});
for j = 1:numel(cases)
 samples(j) = struct('z',cases(j).z0,'u',ctx.uref,'rho',cases(j).rho,'e',cases(j).z0-ctx.zref); %#ok<AGROW>
end
sys = compute_system_matrices(ctx.model,ctx.zref,ctx.uref,ctx.rho_nom,samples);
p = numel(ctx.zref)/2; E = zeros(numel(ctx.zref),p,numel(samples)+1);
E(:,:,1) = ctx.model.fun.B(ctx.zref,ctx.ufull,ctx.rho_nom)*cfg.disturbance_scale;
for j = 1:numel(samples), E(:,:,j+1) = ctx.model.fun.B(samples(j).z,ctx.ufull,samples(j).rho)*cfg.disturbance_scale; end
certificate = certify_design(sys,design.K,design.X,design.Y,design.eta,design.C,design.D,E,cfg,design.kind);
sets = [{'nominal'} {cases.set}];
T = table(sets',certificate.stability_residual-cfg.lmi_margin, ...
 certificate.hinf_residual-cfg.lmi_margin,max(real(certificate.poles),[],1)', ...
 'VariableNames',{'set','unshifted_stability_residual','unshifted_hinf_residual','max_real_pole'});
writetable(T,fullfile(folder,[key '_heldout.csv']));
qsc_atomic_save(fullfile(folder,[key '_heldout.mat']),struct('certificate',certificate,'sys',sys));
end

function run_sampling(cfg,folder)
[ctx,cases,~] = comparison_inputs(cfg); cases = cases(strcmp({cases.scenario},'uncertain_regulation'));
local = cfg; local.active_solver = cfg.comparison_solver; rows = {};
for seed = cfg.seeds
 s = load(fullfile(cfg.output_dir,'comparison',sprintf('QSC_seed%d_design.mat',seed)),'design'); base = s.design;
 for laws = {'uniform_volume','centre_biased','boundary_enriched'}
  law = laws{1}; key = sprintf('%s_seed%d',law,seed); samples = base.samples; n = numel(ctx.zref);
  for j = 1:numel(samples)
   e = base.samples(j).e; rr = norm(e)/cfg.radius; V = min(1,max(0,rr^n));
   if norm(e)==0, continue; end
   radius = rr;
   if strcmp(law,'centre_biased'), radius = V^(4/n);
   elseif strcmp(law,'boundary_enriched')
    if V < 0.5, radius = (2*V)^(1/n); else, radius = (0.9^n+(1-0.9^n)*(2*V-1))^(1/n); end
   end
   samples(j).e = cfg.radius*radius*e/norm(e); samples(j).z = ctx.zref+samples(j).e;
  end
  file = fullfile(folder,[key '_design.mat']);
  if strcmp(law,'uniform_volume')
   design = base; qsc_atomic_save(file,struct('design',design));
  else
   design = qsc_get_design(file,ctx,samples,local,'QSC',seed);
  end
  rows(end+1,:) = {key,'QSC',seed,design.status,design.accepted,design.elapsed}; %#ok<AGROW>
  if design.accepted
   qsc_trials(ctx,design,cases,local,folder,key);
   heldout(ctx,design,cases,local,folder,key);
  end
  write_design_audit(rows,folder);
 end
end
end

function run_saturation(cfg,folder)
[ctx,cases,files] = comparison_inputs(cfg);
for j = 1:numel(files)
 s = load(fullfile(files(j).folder,files(j).name),'design'); d = s.design;
 if ~d.accepted, continue; end
 for limit = cfg.force_limits
  key = sprintf('%s_seed%d_limit%03d',d.kind,d.seed,round(100*limit));
  qsc_trials(ctx,d,cases,cfg,folder,key,limit);
 end
end
end

function run_repeats(stage,cfg,folder)
[ctx,cases,files] = comparison_inputs(cfg);
for j = 1:numel(files)
 s = load(fullfile(files(j).folder,files(j).name),'design'); d = s.design;
 if ~d.accepted, continue; end
 key = sprintf('%s_seed%d',d.kind,d.seed);
 T = qsc_read_trials(fullfile(cfg.output_dir,'comparison',[key '_trials.csv']));
 if strcmp(stage,'force_resolution')
  sub = T(strcmp(T.scenario,'uncertain_regulation') & strcmp(T.set,'uniform'),:);
  [~,worst] = max(sub.peak_commanded_component_N);
  ids = unique([sub.case_id(1);sub.case_id(worst)]);
  for dt = [0.005 0.0005]
   local = cfg; local.dt = dt;
   qsc_trials(ctx,d,cases(ids),local,folder,sprintf('%s_dt%d',key,round(dt*1e6)),Inf,cfg.horizon,true);
  end
 elseif strcmp(stage,'settling')
  ids = select_cases(T); chosen = cases(ids);
  qsc_atomic_save(fullfile(folder,[key '_selection.mat']),struct('case_ids',ids));
  for horizon = [cfg.horizon cfg.long_horizon]
   qsc_trials(ctx,d,chosen,cfg,folder,sprintf('%s_tight%d',key,horizon),Inf,horizon,true);
  end
 else
  assert(exist(fullfile(cfg.output_dir,'saturation','completed.mat'),'file')==2,'Run saturation first');
  for limit = cfg.force_limits
   skey = sprintf('%s_limit%03d',key,round(100*limit));
   T = qsc_read_trials(fullfile(cfg.output_dir,'saturation',[skey '_trials.csv']));
   ids = select_cases(T);
   % Add the most-clipped record of each scenario/set.
   groups = unique(strcat(T.scenario,'|',T.set));
   for g = 1:numel(groups)
    indices = find(strcmp(strcat(T.scenario,'|',T.set),groups{g}));
    [~,worst] = max(T.saturation_duration_s(indices)); ids = [ids;T.case_id(indices(worst))]; %#ok<AGROW>
   end
   ids = unique(ids);
   qsc_atomic_save(fullfile(folder,[skey '_selection.mat']),struct('case_ids',ids));
   qsc_trials(ctx,d,cases(ids),cfg,folder,[skey '_tight100'],limit,cfg.long_horizon,true);
  end
 end
end
end
function ids = select_cases(T)
ids = []; groups = unique(strcat(T.scenario,'|',T.set));
for j = 1:numel(groups)
 sub = T(strcmp(strcat(T.scenario,'|',T.set),groups{j}),:);
 score = sub.final_position_error_m; score(~sub.complete) = Inf;
 [~,worst] = max(score); ids = [ids;sub.case_id(1);sub.case_id(worst)]; %#ok<AGROW>
end
ids = unique(ids);
end

function run_radius_repeats(cfg,folder)
assert(exist(fullfile(cfg.output_dir,'sweeps','completed.mat'),'file')==2,'Run sweeps first');
ctx = qsc_context;
for radius = cfg.radius_values
 sub = fullfile(cfg.output_dir,'sweeps',sprintf('R%03d',round(1000*radius)));
 v = load(fullfile(sub,'validation.mat'),'cases');
 for kind = {'QSC','LQR'}
  seed = cfg.original_design_seed; if strcmp(kind{1},'LQR'), seed = 0; end
  key = sprintf('%s_seed%d',kind{1},seed); s = load(fullfile(sub,[key '_design.mat']),'design'); d = s.design;
  if ~d.accepted, continue; end
  T = qsc_read_trials(fullfile(sub,[key '_trials.csv'])); ids = select_cases(T);
  for horizon = [100 300]
   qsc_trials(ctx,d,v.cases(ids),cfg,folder,sprintf('R%03d_%s_tight%d',round(1000*radius),key,horizon),Inf,horizon,true);
  end
 end
end
end

function run_stabilizability(cfg,folder)
assert(exist(fullfile(cfg.output_dir,'sweeps','completed.mat'),'file')==2,'Run sweeps first');
rows = {}; modes = {}; tols = [1e-10 1e-8 1e-6];
for mask = cfg.masks
 file = fullfile(cfg.output_dir,'sweeps',sprintf('F%02d',mask), ...
  sprintf('QSC_seed%d_design.mat',cfg.original_design_seed));
 s = load(file,'design'); d = s.design;
 A = cat(3,d.sys.A_r,d.sys.A_s); B = cat(3,d.sys.B_r,d.sys.B_s);
 spectra = cell(size(A,3),1);
 for j = 1:size(A,3)
  a = A(:,:,j); b = B(:,:,j); n = size(a,1); C = b; power = b;
  for k = 1:n-1, power = a*power; C = [C power]; end %#ok<AGROW>
  raw = svd(C); norms = sqrt(sum(abs(C).^2,1)); cs = svd(C./max(norms,realmin));
  ranks = arrayfun(@(t)sum(cs>t*cs(1)),tols);
  lambda = eig(a); indices = find(real(lambda)>=-cfg.pole_tolerance);
  passed = true(1,3); pbh = cell(numel(indices),1);
  for k = 1:numel(indices)
   P = [lambda(indices(k))*eye(n)-a b]; sv = svd(P);
   pranks = arrayfun(@(t)sum(sv>t*sv(1)),tols); passed = passed & pranks==n;
   pbh{k} = sv;
   modes(end+1,:) = {mask,j-1,real(lambda(indices(k))),imag(lambda(indices(k))), ...
    sv(end)/sv(1),pranks(1),pranks(2),pranks(3)}; %#ok<AGROW>
  end
  spectra{j} = struct('raw_singular_values',raw,'scaled_singular_values',cs, ...
   'poles',lambda,'PBH_singular_values',{pbh});
  rows(end+1,:) = {mask,j-1,max(real(lambda)),sum(raw>max(size(C))*eps(raw(1))), ...
   ranks(1),ranks(2),ranks(3),numel(indices),passed(1),passed(2),passed(3)}; %#ok<AGROW>
 end
 qsc_atomic_save(fullfile(folder,sprintf('mask%02d_spectra.mat',mask)),struct('spectra',{spectra},'tolerances',tols));
end
T = cell2table(rows,'VariableNames',{'mask','sample_index','max_open_loop_pole', ...
 'raw_rank_default','scaled_rank_1e10','scaled_rank_1e8','scaled_rank_1e6', ...
 'nonstable_modes','PBH_pass_1e10','PBH_pass_1e8','PBH_pass_1e6'});
writetable(T,fullfile(folder,'pair_audit.csv'));
if ~isempty(modes)
 writetable(cell2table(modes,'VariableNames',{'mask','sample_index','real_lambda', ...
  'imag_lambda','PBH_sv_ratio','PBH_rank_1e10','PBH_rank_1e8','PBH_rank_1e6'}),fullfile(folder,'PBH_modes.csv'));
end
end

function run_bounded(cfg,folder)
[ctx,~,~] = comparison_inputs(cfg); rows = {}; cvx_solver(cfg.comparison_solver);
for seed = cfg.seeds
 s = load(fullfile(cfg.output_dir,'comparison',sprintf('QSC_seed%d_design.mat',seed)),'design'); base = s.design;
 sys = base.sys; A = cat(3,sys.A_r,sys.A_s); B = cat(3,sys.B_r,sys.B_s);
 n = size(A,1); m = size(B,2); r0 = cfg.initial_ball_radius; r = cfg.radius;
 for limit = cfg.force_limits
  file = fullfile(folder,sprintf('seed%d_limit%03d.mat',seed,round(100*limit)));
  if exist(file,'file'), saved = load(file,'row'); rows(end+1,:) = saved.row; continue; end %#ok<AGROW>
  h = limit-abs(ctx.uref); started = tic; accepted = false; status = 'NotRun';
  X = []; Y = []; K = []; diagnostics = struct; failure = '';
  try
   assert(all(h>0),'Nonpositive reference-force headroom'); cvx_clear;
   cvx_begin sdp quiet
    variable X(n,n) symmetric
    variable Y(m,n)
    maximize(trace(X))
    X >= cfg.X_floor*eye(n);
    X >= r0^2*eye(n); X <= r^2*eye(n);
    [cfg.Y_bound*eye(m) Y;Y' cfg.Y_bound*eye(n)] >= 0;
    for j = 1:size(A,3)
     H = A(:,:,j)*X+B(:,:,j)*Y; H+H' <= -cfg.lmi_margin*eye(n);
    end
    for j = 1:m, [h(j)^2 Y(j,:);Y(j,:)' X] >= 0; end
   cvx_end
   status = cvx_status;
   if any(strcmp(status,{'Solved','Inaccurate/Solved'}))
    K = Y/X; local = cfg; local.active_solver = cfg.comparison_solver;
    cert = certify_design(sys,K,X,Y,NaN,base.C,base.D,base.E,local,'QSC');
    eigenvalues = eig((X+X')/2); extrema = abs(ctx.uref)+sqrt(max(0,diag(K*X*K')));
    block_eigenvalues = zeros(m,1); extremizer_error = 0;
    for j = 1:m
     block_eigenvalues(j) = min(eig([h(j)^2 Y(j,:);Y(j,:)' X]));
     radius = sqrt(max(0,K(j,:)*X*K(j,:)'));
     if radius>0
      extremizer = X*K(j,:)'/radius;
      explicit = max(abs([ctx.uref(j)+K(j,:)*extremizer,ctx.uref(j)-K(j,:)*extremizer]));
      extremizer_error = max(extremizer_error,abs(explicit-extrema(j)));
     end
    end
    diagnostics = struct('certificate',cert,'initial_ball_residual',max(eig(r0^2*eye(n)-X)), ...
     'support_residual',max(eig(X-r^2*eye(n))),'bound_min_eigenvalues',block_eigenvalues, ...
     'command_maxima',extrema,'axes',sqrt(eigenvalues), ...
     'volume_ratio',exp(0.5*sum(log(eigenvalues))-n*log(r)), ...
     'extremizer_error_N',extremizer_error);
    accepted = cert.accepted && diagnostics.initial_ball_residual<=cfg.X_tolerance ...
     && diagnostics.support_residual<=cfg.X_tolerance ...
     && min(block_eigenvalues)>=-cfg.residual_tolerance && max(extrema-limit)<=1e-7;
   end
  catch exception
   status = 'Exception'; failure = getReport(exception,'basic','hyperlinks','off');
   cvx_clear;
  end
  elapsed = toc(started); amin = NaN; amax = NaN; vol = NaN; peak = NaN; residual = NaN;
  if isfield(diagnostics,'axes')
   amin = min(diagnostics.axes); amax = max(diagnostics.axes); vol = diagnostics.volume_ratio;
   peak = max(diagnostics.command_maxima); residual = diagnostics.certificate.max_stability_residual;
  end
  row = {seed,limit,status,accepted,elapsed,min(h),amin,amax,vol,peak,residual,failure};
  qsc_atomic_save(file,struct('X',X,'Y',Y,'K',K,'diagnostics',diagnostics, ...
   'samples',base.samples,'sys',sys,'reference',rmfield(ctx,'model'),'row',{row}));
  rows(end+1,:) = row; %#ok<AGROW>
  fprintf('BOUNDED seed%d limit%.2f: %s accepted=%d\n',seed,limit,status,accepted);
 end
end
writetable(cell2table(rows,'VariableNames',{'seed','limit_N','status','accepted', ...
 'synthesis_time_s','minimum_headroom_N','minimum_axis','maximum_axis','volume_ratio', ...
 'peak_command_bound_N','training_residual','failure'}),fullfile(folder,'controller_audit.csv'));
end

function run_margins(stage,cfg,folder)
rows = {}; local = cfg; local.active_solver = cfg.comparison_solver;
if strcmp(stage,'margin')
 grids = 2; counts = cfg.N; scales = cfg.margin_scales; validation_count = cfg.margin_count;
else
 grids = cfg.scaling_grids; counts = cfg.scaling_N; scales = 1; validation_count = cfg.scaling_validation_count;
end
for grid = grids
 ctx = qsc_context(15,grid);
 for N = counts
  samples = generate_samples(ctx.zref,ctx.uref,ctx.rho_nom,cfg.uncertainty,N,cfg.radius,2,1);
  key = sprintf('grid%d_N%d',grid,N);
  design = qsc_get_design(fullfile(folder,[key '_design.mat']),ctx,samples,local,'QSC',2);
  for scale = scales
   file = fullfile(folder,sprintf('%s_scale%03d.mat',key,round(100*scale)));
   if exist(file,'file')
    s = load(file,'row'); rows(end+1,:) = s.row; continue; %#ok<AGROW>
   end
   if strcmp(stage,'margin'), seed = 12346+(scale>1); else, seed = 54322+grid-1; end
   validation = generate_samples(ctx.zref,ctx.uref,ctx.rho_nom,cfg.uncertainty, ...
    validation_count,scale*cfg.radius,seed,1);
   margins = nan(validation_count,1); distances = zeros(validation_count,1);
   if design.accepted
    P = design.X\eye(numel(ctx.zref));
    for j = 1:validation_count
     v = validation(j); distances(j) = norm(v.e);
     Ac = ctx.model.fun.A(v.z,ctx.uref,v.rho)+ctx.model.fun.B_mod(v.z,ctx.uref,v.rho)*design.K;
     M = Ac'*P+P*Ac; margins(j) = max(eig((M+M')/2));
    end
   end
   pass_fraction = NaN;
   if design.accepted, pass_fraction = mean(margins<0); end
   row = {grid,numel(ctx.zref),N,scale,design.status,design.accepted,design.elapsed, ...
    validation_count,sum(margins<0),pass_fraction,max(margins)};
   qsc_atomic_save(file,struct('row',{row},'samples',validation,'margins',margins));
   writetable(table((1:validation_count)',distances,margins, ...
    'VariableNames',{'case_id','radial_distance','lyapunov_margin'}),strrep(file,'.mat','.csv'));
   rows(end+1,:) = row; %#ok<AGROW>
  end
 end
end
writetable(cell2table(rows,'VariableNames',{'grid','state_dimension','N','radius_scale', ...
 'status','accepted','synthesis_time_s','validation_count','certificate_passes', ...
 'heldout_pass_fraction','worst_margin'}),fullfile(folder,'summary.csv'));
end

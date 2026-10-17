function design = qsc_get_design(file,ctx,samples,cfg,kind,seed)
if exist(file,'file'), data = load(file,'design'); design = data.design; return; end
cvx_clear;
cvx_solver(cfg.active_solver);
design = synthesize_design(ctx,samples,cfg,kind);
design.seed = seed;
design.reference = rmfield(ctx,'model');
qsc_atomic_save(file,struct('design',design));
fprintf('DESIGN %s seed%d: %s accepted=%d time=%.3fs\n', ...
 kind,seed,design.status,design.accepted,design.elapsed);
end

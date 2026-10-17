function cfg = qsc_settings(mode)
% Explicit, editable protocol. Quick mode is an installation check only.
if nargin < 1, mode = 'quick'; end
assert(any(strcmp(mode, {'quick','full'})), 'Use quick or full');
cfg.schema = 1;
cfg.mode = mode;
cfg.cvx_startup = ''; % Optional full path to your cvx_startup.m.
cfg.comparison_solver = 'sdpt3';
cfg.original_solver = 'sedumi';
cfg.N = 100; cfg.radius = 0.05; cfg.uncertainty = 0.1;
cfg.seeds = 2:6; cfg.validation_seed = 22001;
cfg.uniform_count = 100; cfg.stress_count = 32;
cfg.original_validation_count = 1000;
cfg.original_validation_seed = 101;
cfg.original_design_seed = 2;
cfg.N_values = [50 75 100 125 150 200];
cfg.radius_values = [0.1 0.075 0.05 0.025];
cfg.uncertainty_values = [0.2 0.15 0.1 0.05];
cfg.diagnostic_uncertainties = [0.05 0.1 0.15];
cfg.masks = 1:15;
cfg.position_scale = 0.05; cfg.force_scale = 1;
cfg.disturbance_scale = 1;
cfg.lmi_margin = 1e-4; cfg.X_floor = 1e-6; cfg.Y_bound = 8;
cfg.residual_tolerance = 1e-7; cfg.X_tolerance = 1e-8;
cfg.pole_tolerance = 1e-8;
cfg.horizon = 30; cfg.long_horizon = 100;
cfg.original_horizon = 100; cfg.event_threshold = 1e-8;
cfg.dt = 0.002; cfg.reltol = 1e-7; cfg.abstol = 1e-9;
cfg.maxstep = 0.05; cfg.tight_reltol = 1e-9; cfg.tight_abstol = 1e-11;
cfg.settle_velocity = 1e-6; cfg.settle_acceleration = 1e-6;
cfg.settle_position_range = 1e-5;
cfg.force_limits = [0.3 0.75 2];
cfg.initial_ball_radius = 0.01;
cfg.margin_count = 5000; cfg.margin_scales = [1 1.5];
cfg.scaling_grids = 2:4; cfg.scaling_N = [25 50 75 100 150 200 300];
cfg.scaling_validation_count = 1000;
cfg.block_size = 20;
cfg.save_traces = false; % Metrics use ALL integration points regardless.
cfg.trace_stride = 25;
cfg.stages = {'mapping','sweeps','diagnostic','comparison','sampling', ...
 'saturation','settling','force_resolution','radius_settling', ...
 'saturation_settling','stabilizability','bounded','margin','scaling'};
if strcmp(mode, 'quick')
 cfg.N = 5; cfg.seeds = 2;
 cfg.uniform_count = 2; cfg.stress_count = 2;
 cfg.original_validation_count = 2;
 cfg.N_values = 5; cfg.radius_values = 0.05;
 cfg.uncertainty_values = 0.1; cfg.diagnostic_uncertainties = 0.05;
 cfg.masks = [5 15]; cfg.mapping_masks = [1 5 15];
 cfg.margin_count = 20; cfg.scaling_grids = 2;
 cfg.scaling_N = 5; cfg.scaling_validation_count = 20;
 cfg.dt = 0.02;
 cfg.stages = {'mapping','comparison'};
else
 cfg.mapping_masks = 1:15;
end
end

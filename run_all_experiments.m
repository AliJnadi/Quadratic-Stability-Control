function run_all_experiments(mode, overrides)
% Quick mode resets results/quick; full mode resumes matching saved runs.
% run_all_experiments('quick')
% run_all_experiments('full', struct('cvx_startup', '/path/cvx_startup.m'))
if nargin < 1, mode = 'quick'; end
if nargin < 2, overrides = struct; end
root = fileparts(mfilename('fullpath'));
previous_path = path; cleanup = onCleanup(@()path(previous_path)); %#ok<NASGU>
addpath(fullfile(root,'src'), fullfile(root,'experiments'), ...
 fullfile(root,'analysis'), '-begin');
cfg = qsc_settings(mode);
names = fieldnames(overrides);
for j = 1:numel(names)
 assert(isfield(cfg,names{j}) || strcmp(names{j},'output_dir'), ...
  'Unknown setting: %s', names{j});
 cfg.(names{j}) = overrides.(names{j});
end
if ~isfield(cfg,'output_dir'), cfg.output_dir = fullfile(root,'results',mode); end
output_file = java.io.File(cfg.output_dir);
cfg.output_dir = char(output_file.getCanonicalPath());
assert(cfg.block_size > 0 && cfg.dt > 0 && cfg.N >= 1);
assert(cfg.X_floor == 1e-6 && cfg.Y_bound == 8, ...
 'The retained QSC solver uses X_floor=1e-6 and Y_bound=8.');
assert(~isempty(cfg.seeds) && ~isempty(cfg.stages));
if ~isempty(cfg.cvx_startup), run(cfg.cvx_startup); end
assert(exist('cvx_begin','file') == 2, 'Install CVX or set cfg.cvx_startup');
assert(exist('sym','file') == 2, 'Symbolic Math Toolbox is required');
assert(exist('lqr','file') == 2, 'Control System Toolbox is required');
% CVX startup can alter paths; explicitly restore project precedence.
addpath(fullfile(root,'src'),'-begin');
clear qsc_context; % Discard stale symbolic handles from earlier calls.
assert(strcmp(which('expand_input'),fullfile(root,'src','expand_input.m')));
cfg.active_solver = cfg.comparison_solver;
% Reset only the standard quick output, after dependency checks have passed.
% Custom output folders and all full runs retain their resume behavior.
root_file = java.io.File(root);
standard_quick = fullfile(char(root_file.getCanonicalPath()),'results','quick');
if strcmp(mode,'quick') && strcmp(cfg.output_dir,standard_quick)
 if exist(cfg.output_dir,'dir')
  fprintf('QUICK: clearing previous results in %s\n',cfg.output_dir);
  [removed,message] = rmdir(cfg.output_dir,'s');
  assert(removed,'Could not clear quick results: %s',message);
 end
end
if ~exist(cfg.output_dir,'dir'), mkdir(cfg.output_dir); end
manifest = qsc_manifest(root,cfg);
file = fullfile(cfg.output_dir,'manifest.mat');
if exist(file,'file')
 old = load(file,'manifest');
 assert(strcmp(old.manifest.signature,manifest.signature), ...
  'Code/settings/environment changed. Choose a new output_dir; do not mix runs.');
else
 qsc_atomic_save(file,struct('manifest',manifest,'cfg',cfg));
 qsc_write_json(fullfile(cfg.output_dir,'manifest.json'),manifest);
end
diary(fullfile(cfg.output_dir,'run.log'));
log_cleanup = onCleanup(@()diary('off')); %#ok<NASGU>
fprintf('QSC %s mode: %d stages. Results folder: %s\n',mode,numel(cfg.stages),cfg.output_dir);
for j = 1:numel(cfg.stages)
 stage = cfg.stages{j}; folder = fullfile(cfg.output_dir,stage);
 if ~exist(folder,'dir'), mkdir(folder); end
 flag = fullfile(folder,'completed.mat');
 if exist(flag,'file'), fprintf('RESUME: %s already complete\n',stage); continue; end
 fprintf('\nSTART %s (%d/%d) %s\n',stage,j,numel(cfg.stages),datestr(now,30));
 started = tic;
 try
  qsc_run_stage(stage,cfg,folder);
  qsc_atomic_save(flag,struct('stage',stage,'elapsed_s',toc(started), ...
   'finished',datestr(now,30),'signature',manifest.signature));
 catch exception
  qsc_write_json(fullfile(folder,'failure.json'),struct('stage',stage, ...
   'message',getReport(exception,'extended','hyperlinks','off')));
  rethrow(exception); % Keep completed blocks; never mark an interrupted stage complete.
 end
end
qsc_generate_outputs(cfg.output_dir);
fprintf('\nFinished. Inspect tables, figures and acceptance/failure counts before use.\n');
end

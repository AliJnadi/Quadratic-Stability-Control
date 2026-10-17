function qsc_generate_outputs(output_dir)
% Regenerate outputs without re-running any experiment.
output_file = java.io.File(output_dir);
output_dir = char(output_file.getCanonicalPath());
root = fileparts(fileparts(mfilename('fullpath')));
old_path = path; cleanup = onCleanup(@()path(old_path)); %#ok<NASGU>
addpath(fullfile(root,'src'));
folders = dir(output_dir); status_rows = {};
for f = 1:numel(folders)
 if ~folders(f).isdir || folders(f).name(1)=='.' || strcmp(folders(f).name,'verification'), continue; end
 stage = folders(f).name; folder = fullfile(output_dir,stage);
 files = dir(fullfile(folder,'**','*_trials.csv'));
 rows = {};
 for j = 1:numel(files)
  file = fullfile(files(j).folder,files(j).name); T = qsc_read_trials(file);
  [~,key] = fileparts(file); key = key(1:end-7); % Remove _trials.
  experiment = regexprep(key,'_seed[0-9]+','');
  if any(strcmp(stage,{'sweeps','diagnostic'}))
   [~,setting] = fileparts(files(j).folder);
   experiment = [setting '_' experiment];
   key = [setting '_' key];
  end
  groups = unique(strcat(T.scenario,'|',T.set));
  for g = 1:numel(groups)
   sub = T(strcmp(strcat(T.scenario,'|',T.set),groups{g}),:);
   ok = logical(sub.complete); settled = NaN;
   if ismember('settled',sub.Properties.VariableNames), settled = sum(sub.settled); end
   limit = Inf;
   if ismember('limit_N',sub.Properties.VariableNames), limit = sub.limit_N(1); end
   rows(end+1,:) = {experiment,key,sub.controller{1},sub.seed(1),sub.scenario{1}, ...
    sub.set{1},limit,height(sub),sum(ok),settled, ...
    mean(sub.position_MSE_m2(ok)),std(sub.position_MSE_m2(ok)), ...
    mean(sub.final_position_error_m(ok)),std(sub.final_position_error_m(ok)), ...
    mean(sub.incremental_force_RMS_N(ok)),std(sub.incremental_force_RMS_N(ok)), ...
    max_or_nan(sub.peak_commanded_component_N(ok))}; %#ok<AGROW>
  end
 end
 if ~isempty(rows)
  seeds = cell2table(rows,'VariableNames',{'experiment','key','controller','seed','scenario', ...
   'set','limit_N','trials','completed','settled','position_MSE_mean','position_MSE_case_SD', ...
   'final_error_mean','final_error_case_SD','force_RMS_mean','force_RMS_case_SD','peak_command_N'});
  writetable(seeds,fullfile(folder,'seed_summary.csv'));
  summary = aggregate_seeds(seeds);
  writetable(summary,fullfile(folder,'aggregate_summary.csv'));
  write_latex(summary,fullfile(folder,'results_table.tex'));
  status_rows(end+1,:) = {stage,sum(seeds.trials),sum(seeds.completed), ...
   sum(seeds.trials-seeds.completed)}; %#ok<AGROW>
  if any(strcmp(stage,{'comparison','sampling','saturation'})), tradeoff_plot(summary,folder,stage); end
  if any(strcmp(stage,{'comparison','sampling','saturation'})), detail_plot(files,folder,stage); end
  if any(strcmp(stage,{'settling','force_resolution','radius_settling','saturation_settling'}))
   repeat_report(files,folder,stage,output_dir);
  end
 end
 if any(strcmp(stage,{'sweeps','diagnostic'})) && exist(fullfile(folder,'summary.csv'),'file')
  qsc_export_legacy_summary(folder);
  sweep_plot(folder,stage);
 elseif strcmp(stage,'scaling') && exist(fullfile(folder,'summary.csv'),'file')
  scaling_plot(folder);
 elseif strcmp(stage,'margin') && exist(fullfile(folder,'summary.csv'),'file')
  margin_plot(folder);
 elseif strcmp(stage,'stabilizability') && exist(fullfile(folder,'pair_audit.csv'),'file')
  ranks_plot(folder);
 elseif strcmp(stage,'bounded') && exist(fullfile(folder,'controller_audit.csv'),'file')
  bounded_plot(folder);
 end
end
if ~isempty(status_rows)
 T = cell2table(status_rows,'VariableNames',{'stage','trial_records','completed','incomplete'});
 writetable(T,fullfile(output_dir,'RUN_COUNTS.csv'));
end
fid = fopen(fullfile(output_dir,'RESULTS.md'),'w'); assert(fid>=0);
fprintf(fid,'# Locally generated QSC results\n\nGenerated %s.\n\n',datestr(now,30));
fprintf(fid,'Inspect RUN_COUNTS.csv, controller_audit.csv and all per-trial failures before using results. A completed stage does not mean every design was accepted or every trajectory settled.\n\n');
fprintf(fid,'seed_summary.csv contains case means for each gain; aggregate_summary.csv averages these gain-level means. Between-seed SD is not a confidence interval. Repeated stress tests and settling repeats are not independent evidence.\n\n');
fprintf(fid,'CSV trial records and MAT designs remain local. Figures and tables are generated from these results, not from historical published values. Hardware timings and feasible QSC gains may change with solver versions.\n\n');
fprintf(fid,'The full continuous nonlinear domain is not certified. Saturation tests are illustrative; bounded-command results do not prove nonlinear invariance. Uncertain disturbance ratios include reference drift.\n');
fclose(fid);
end
function result = max_or_nan(v)
if isempty(v), result = NaN; else, result = max(v); end
end
function summary = aggregate_seeds(T)
keys = strcat(T.experiment,'|',T.scenario,'|',T.set);
groups = unique(keys); rows = {};
for j = 1:numel(groups)
 sub = T(strcmp(keys,groups{j}),:);
 rows(end+1,:) = {sub.experiment{1},sub.controller{1},sub.scenario{1},sub.set{1}, ...
  sub.limit_N(1),height(sub),sum(sub.trials),sum(sub.completed), ...
  mean(sub.position_MSE_mean),std(sub.position_MSE_mean), ...
  mean(sub.final_error_mean),std(sub.final_error_mean), ...
  mean(sub.force_RMS_mean),std(sub.force_RMS_mean),max(sub.peak_command_N)}; %#ok<AGROW>
end
summary = cell2table(rows,'VariableNames',{'experiment','controller','scenario','set','limit_N', ...
 'gains','trials','completed','position_MSE_mean','position_MSE_seed_SD', ...
 'final_error_mean','final_error_seed_SD','force_RMS_mean','force_RMS_seed_SD','peak_command_N'});
end
function write_latex(T,file)
fid = fopen(file,'w'); assert(fid>=0); cleanup = onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%% Automatically generated; verify completion counts before manuscript use.\n');
fprintf(fid,'\\begin{tabular}{llllrrr}\n\\hline\nDesign & Scenario & Set & Gains & Position MSE & Final error & Force RMS \\\\\n\\hline\n');
for j = 1:height(T)
 fprintf(fid,'%s & %s & %s & %d & %.6g & %.6g & %.6g \\\\\n', ...
  escape_tex(T.experiment{j}),escape_tex(T.scenario{j}),escape_tex(T.set{j}), ...
  T.gains(j),T.position_MSE_mean(j),T.final_error_mean(j),T.force_RMS_mean(j));
end
fprintf(fid,'\\hline\n\\end{tabular}\n');
end
function t = escape_tex(t)
t = strrep(t,'_','\_');
end
function finish_figure(fig,folder,name)
set(fig,'Color','w','Units','inches');
figure_position = get(fig,'Position');
page_size = figure_position(3:4);
set(fig,'PaperUnits','inches','PaperSize',page_size, ...
 'PaperPosition',[0 0 page_size],'PaperPositionMode','manual');
print(fig,fullfile(folder,[name '.pdf']),'-dpdf','-painters');
print(fig,fullfile(folder,[name '.png']),'-dpng','-r160'); close(fig);
end
function tradeoff_plot(T,folder,stage)
sub = T(strcmp(T.scenario,'uncertain_regulation'),:);
fig = figure('Visible','off','Position',[100 100 850 430]); hold on;
for j = 1:height(sub)
 if sub.completed(j) ~= sub.trials(j), continue; end
 label = [sub.experiment{j} ' / ' sub.set{j}];
 loglog(sub.force_RMS_mean(j),sub.position_MSE_mean(j),'o','MarkerSize',7,'LineWidth',1.4,'DisplayName',label);
end
set(gca,'XScale','log','YScale','log'); grid on;
xlabel('Incremental requested force RMS [N]'); ylabel('Position time-weighted MSE [m^2]');
legend('show','Interpreter','none','Location','eastoutside'); title([stage ' error-effort comparison'],'Interpreter','none');
finish_figure(fig,folder,'error_effort');
end
function sweep_plot(folder,stage)
T = readtable(fullfile(folder,'summary.csv'));
T = T(strcmp(T.scenario,'uncertain_regulation'),:);
for sweep = unique(T.sweep)'
 S = T(strcmp(T.sweep,sweep{1}),:); fig = figure('Visible','off','Position',[100 100 1000 330]);
 labels = {'Position MSE [m^2]','Final position error [m]','Incremental force RMS [N]'};
 columns = {'position_MSE_mean','final_position_error_mean','incremental_force_RMS_mean'};
 for panel = 1:3
  subplot(1,3,panel); hold on;
  for kind = {'QSC','LQR'}
   K = S(strcmp(S.controller,kind{1}) & S.completed==S.trials,:);
   values = unique(K.value); means = zeros(size(values)); spreads = means;
   for j = 1:numel(values)
    v = K.(columns{panel})(K.value==values(j)); means(j) = mean(v); spreads(j) = std(v);
   end
   errorbar(values,means,spreads,'-o','DisplayName',kind{1});
  end
  xlabel(sweep{1}); ylabel(labels{panel}); grid on; legend('show','Location','best');
  title([stage ' / ' sweep{1}],'Interpreter','none');
 end
 finish_figure(fig,folder,['sweep_' sweep{1}]);
end
end
function scaling_plot(folder)
T = readtable(fullfile(folder,'summary.csv')); fig = figure('Visible','off','Position',[100 100 850 350]);
for panel = 1:2
 subplot(1,2,panel); hold on;
 for grid_size = unique(T.grid)'
  sub = T(T.grid==grid_size,:); [~,order] = sort(sub.N); sub = sub(order,:);
  if panel==1, y = sub.synthesis_time_s; else, y = sub.heldout_pass_fraction; end
  plot(sub.N,y,'-o','DisplayName',sprintf('n=%d',sub.state_dimension(1)));
 end
 xlabel('Design samples N'); grid on; legend('show','Location','best');
 if panel==1, ylabel('Observed synthesis time [s]'); else, ylabel('Held-out certificate pass fraction'); end
end
finish_figure(fig,folder,'scaling');
end
function margin_plot(folder)
files = dir(fullfile(folder,'grid*_scale*.csv')); fig = figure('Visible','off'); hold on;
for j = 1:numel(files)
 T = readtable(fullfile(files(j).folder,files(j).name));
 plot(T.radial_distance,T.lyapunov_margin,'.','DisplayName',files(j).name);
end
xlabel('Numerical distance from reference'); ylabel('Largest Lyapunov-derivative eigenvalue');
grid on; legend('show','Interpreter','none','Location','best');
finish_figure(fig,folder,'heldout_margins');
end
function ranks_plot(folder)
T = readtable(fullfile(folder,'pair_audit.csv')); masks = unique(T.mask); values = zeros(numel(masks),2);
for j = 1:numel(masks)
 sub = T(T.mask==masks(j),:); values(j,:) = [min(sub.scaled_rank_1e8) max(sub.scaled_rank_1e8)];
end
fig = figure('Visible','off'); bar(masks,values); xlabel('Actuator mask'); ylabel('Numerical rank');
legend('Minimum','Maximum','Location','best'); grid on;
finish_figure(fig,folder,'stabilizability_ranks');
end
function bounded_plot(folder)
T = readtable(fullfile(folder,'controller_audit.csv')); fig = figure('Visible','off','Position',[100 100 850 350]);
subplot(1,2,1); plot(T.limit_N,T.volume_ratio,'o'); xlabel('Illustrative component limit [N]'); ylabel('Ellipsoid/design-ball volume ratio'); grid on;
subplot(1,2,2); plot(T.limit_N,T.peak_command_bound_N,'o'); hold on;
plot(T.limit_N,T.limit_N,'k--'); xlabel('Illustrative component limit [N]'); ylabel('Verified command bound [N]'); grid on;
finish_figure(fig,folder,'bounded_command_regions');
end
function detail_plot(files,folder,stage)
rows = {};
for j = 1:numel(files)
 T = qsc_read_trials(fullfile(files(j).folder,files(j).name));
 [~,key] = fileparts(files(j).name); key = key(1:end-7);
 if strcmp(stage,'comparison')
  sub = T(strcmp(T.scenario,'nominal_disturbance'),:);
  if isempty(sub), continue; end
  if any(~sub.complete), continue; end
  design = load(fullfile(folder,[key '_design.mat']),'design');
  rows(end+1,:) = {key,mean(sub.weighted_L2_ratio),design.design.elapsed,norm(design.design.K,'fro')}; %#ok<AGROW>
 else
  sub = T(strcmp(T.scenario,'uncertain_regulation') & strcmp(T.set,'uniform'),:);
  if isempty(sub), continue; end
  if any(~sub.complete), continue; end
  clipping = 0;
  if ismember('saturation_duration_s',sub.Properties.VariableNames), clipping = mean(sub.saturation_duration_s); end
  rows(end+1,:) = {key,mean(sub.position_MSE_m2),mean(sub.incremental_force_RMS_N),clipping}; %#ok<AGROW>
 end
end
if isempty(rows), return; end
T = cell2table(rows,'VariableNames',{'key','metric1','metric2','metric3'});
writetable(T,fullfile(folder,'figure_metrics.csv'));
fig = figure('Visible','off','Position',[100 100 1100 470]);
labels = {'Position MSE [m^2]','Incremental requested force RMS [N]','Clipping duration [s]'};
if strcmp(stage,'comparison'), labels = {'Nominal empirical normalized L2 ratio','Observed synthesis time [s]','Gain Frobenius norm'}; end
for panel = 1:3
 subplot(1,3,panel); bar(T.(sprintf('metric%d',panel)));
 set(gca,'XTick',1:height(T),'XTickLabel',T.key,'TickLabelInterpreter','none');
 xtickangle(65); ylabel(labels{panel}); grid on;
end
finish_figure(fig,folder,'performance_details');
end

function repeat_report(files,folder,stage,output_dir)
rows = {};
for j = 1:numel(files)
 T = qsc_read_trials(fullfile(files(j).folder,files(j).name));
 [~,key] = fileparts(files(j).name); key = key(1:end-7);
 basekey = sprintf('%s_seed%d',T.controller{1},T.seed(1));
 baseline_folder = fullfile(output_dir,'comparison');
 if strcmp(stage,'radius_settling')
  setting = regexp(key,'^R[0-9]+','match','once');
  baseline_folder = fullfile(output_dir,'sweeps',setting);
 elseif strcmp(stage,'saturation_settling')
  basekey = regexp(key,'^.*_seed[0-9]+_limit[0-9]+','match','once');
  baseline_folder = fullfile(output_dir,'saturation');
 end
 baseline = qsc_read_trials(fullfile(baseline_folder,[basekey '_trials.csv']));
 for k = 1:height(T)
  b = baseline(baseline.case_id==T.case_id(k),:);
  assert(height(b)==1,'Repeat case must have one baseline record');
  rows(end+1,:) = {key,T.case_id(k),logical(b.complete),logical(T.complete(k)), ...
   abs(T.final_position_error_m(k)-b.final_position_error_m), ...
   abs(T.position_MSE_m2(k)-b.position_MSE_m2), ...
   100*abs(T.incremental_force_RMS_N(k)-b.incremental_force_RMS_N)/max(abs(b.incremental_force_RMS_N),eps), ...
   logical(T.settled(k))}; %#ok<AGROW>
 end
end
if isempty(rows), return; end
T = cell2table(rows,'VariableNames',{'key','case_id','baseline_complete','repeat_complete', ...
 'absolute_final_error_change_m','absolute_position_MSE_change_m2','force_RMS_relative_change_percent','repeat_settled'});
writetable(T,fullfile(folder,'repeat_differences.csv'));
fig = figure('Visible','off','Position',[100 100 900 350]);
subplot(1,2,1); plot(T.absolute_final_error_change_m,'.');
xlabel('Selected repeat record'); ylabel('Absolute final-error change [m]'); grid on;
subplot(1,2,2); plot(T.force_RMS_relative_change_percent,'.');
xlabel('Selected repeat record'); ylabel('Force RMS relative change [%]'); grid on;
finish_figure(fig,folder,'repeat_checks');
end

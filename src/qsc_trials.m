function T = qsc_trials(ctx,design,cases,cfg,folder,key,limit,horizon,tight,legacy)
if nargin < 7, limit = Inf; end
if nargin < 8, horizon = cfg.horizon; end
if nargin < 9, tight = false; end
if nargin < 10, legacy = false; end
if ~exist(folder,'dir'), mkdir(folder); end
protocol_file = fullfile(folder,[key '_protocol.mat']);
if ~exist(protocol_file,'file')
 qsc_atomic_save(protocol_file,struct('cases',cases,'cfg',cfg,'key',key, ...
  'limit',limit,'horizon',horizon,'tight',tight,'legacy',legacy));
end
T = table;
for first = 1:cfg.block_size:numel(cases)
 last = min(first+cfg.block_size-1,numel(cases));
 file = fullfile(folder,sprintf('%s_block%05d.mat',key,first));
 if exist(file,'file')
  data = load(file,'metrics'); metrics = data.metrics;
 else
  metrics = []; traces = cell(0,1);
  for j = first:last
   if legacy
    m = qsc_legacy_trial(ctx,design,cases(j),cfg,horizon);
   else
    r = simulate_case(ctx,design,cases(j),cfg,limit,horizon,tight,false);
    m = r.metrics; m.seed = design.seed;
    if cfg.save_traces
     traces{end+1,1} = struct('case_id',cases(j).id,'trace',r.trace); %#ok<AGROW>
    end
   end
   metrics = [metrics; m]; %#ok<AGROW>
  end
  qsc_atomic_save(file,struct('metrics',metrics,'traces',{traces}));
 end
 block = struct2table(metrics);
 csv_file = strrep(file,'.mat','.csv');
 if ~exist(csv_file,'file'), writetable(block,csv_file); end
 T = [T; block]; %#ok<AGROW>
 fprintf('%s %d/%d trials; block failures=%d\n',key,last,numel(cases),sum(~block.complete));
end
writetable(T,fullfile(folder,[key '_trials.csv']));
end

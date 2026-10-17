function summary = qsc_export_legacy_summary(folder)
% Like-for-like metrics for the original Final_Table; retain physical outputs.
files = dir(fullfile(folder,'*','*_trials.csv')); rows = {};
for j = 1:numel(files)
 T = qsc_read_trials(fullfile(files(j).folder,files(j).name));
 if ~ismember('legacy_mixed_state_MSE',T.Properties.VariableNames), continue; end
 [~,setting] = fileparts(files(j).folder);
 for scenario = {'nominal_regulation','uncertain_regulation'}
  S = T(strcmp(T.scenario,scenario{1}),:); if isempty(S), continue; end
  ok = logical(S.complete);
  rows(end+1,:) = {setting,S.controller{1},S.seed(1),scenario{1},height(S),sum(ok), ...
   mean(S.legacy_mixed_state_MSE(ok)),std(S.legacy_mixed_state_MSE(ok)), ...
   mean(S.legacy_mixed_final_norm(ok)),std(S.legacy_mixed_final_norm(ok))}; %#ok<AGROW>
 end
end
summary = cell2table(rows,'VariableNames',{'setting','controller','seed','scenario', ...
 'trials','completed','legacy_mixed_MSE_mean','legacy_mixed_MSE_case_SD', ...
 'legacy_mixed_final_norm_mean','legacy_mixed_final_norm_case_SD'});
if ~isempty(rows), writetable(summary,fullfile(folder,'legacy_summary.csv')); end
end

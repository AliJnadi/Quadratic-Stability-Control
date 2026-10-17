function qsc_case_provenance(cases,folder)
% Identify repeated physical inputs; retain original case weights and IDs.
first = zeros(numel(cases),1); keys = cell(numel(cases),1);
for j = 1:numel(cases)
 keys{j} = jsonencode(struct('initial_state',cases(j).z0, ...
  'parameters',cases(j).rho,'disturbed',cases(j).disturbed));
 previous = find(strcmp(keys(1:j-1),keys{j}),1);
 first(j) = j; if ~isempty(previous), first(j) = previous; end
end
T = table([cases.id]',{cases.set}',{cases.scenario}',first, ...
 first~=(1:numel(cases))','VariableNames', ...
 {'case_id','set','scenario','first_identical_input_case','duplicate_input'});
writetable(T,fullfile(folder,'case_provenance.csv'));
end

function T = qsc_read_trials(file)
options = detectImportOptions(file);
for name = {'failure','controller','set','scenario'}
 if ismember(name{1},options.VariableNames)
  options = setvartype(options,name{1},'char');
 end
end
T = readtable(file,options);
end

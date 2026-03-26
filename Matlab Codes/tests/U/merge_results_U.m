% Folder containing the .mat files
files = dir(fullfile('*.mat'));

col_name = 'SamplingRadius';

Tall = table();  % final combined table

for k = 1:numel(files)

    % Load file
    filePath = fullfile(files(k).folder, files(k).name);
    S = load(filePath);

    % Extract table T
    if ~isfield(S, 'T')
        warning('No table T in %s, skipping.', files(k).name);
        continue
    end
    T = S.T;

    % -------------------------------------------------
    % Extract number of samples from file name or variable name
    % Example: N100R_05U10.mat  --> 100
    % Uncomment the one which is related to your data       
    % -------------------------------------------------
    tokens = regexp(files(k).name, 'N100R50U(\d+)', 'tokens');
    
    Nsamples = str2double(tokens{1}{1})*1e-1;

    % -------------------------------------------------
    % Add number of samples column (after Method)
    % -------------------------------------------------
    Ncol = repmat(Nsamples, height(T), 1);
    T = addvars(T, Ncol, 'After', 'Method', ...
        'NewVariableNames', col_name);

    % Append to master table
    Tall = [Tall; T];
end
Tall = sortrows(Tall, col_name);
% Display final table
disp(Tall)
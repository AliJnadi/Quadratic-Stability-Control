% % Folder containing the .mat files
% files = dir(fullfile('*.mat'));
% 
% col_name = 'NumSamples';
% 
% Tall = table();  % final combined table
% 
% for k = 1:numel(files)
% 
%     % Load file
%     filePath = fullfile(files(k).folder, files(k).name);
%     S = load(filePath);
% 
%     % Extract table T
%     if ~isfield(S, 'T')
%         warning('No table T in %s, skipping.', files(k).name);
%         continue
%     end
%     T = S.T;
%     
%     % -------------------------------------------------
%     % Extract number of samples from file name or variable name
%     % Example: N100R_05U10.mat  --> 100
%     % Uncomment the one which is related to your data       
%     % -------------------------------------------------
%     tokens = regexp(files(k).name, 'N(\d+)R', 'tokens');
%     
%     if isempty(tokens)
%         warning('Cannot extract N from %s, skipping.', files(k).name);
%         continue
%     end
%     
%     Nsamples = str2double(tokens{1}{1});
% 
%     % -------------------------------------------------
%     % Add number of samples column (after Method)
%     % -------------------------------------------------
%     Ncol = repmat(Nsamples, height(T), 1);
%     T = addvars(T, Ncol, 'After', 'Method', ...
%         'NewVariableNames', col_name);
%     
%     
%     % Append to master table
%     Tall = [Tall; T];
% end
% Tall = sortrows(Tall, col_name);
% % Display final table
% disp(Tall)

% Folder containing the .mat files
files = dir(fullfile('*.mat'));

col_name = 'NumSamples';
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
    
    % Extract number of samples from filename
    tokens = regexp(files(k).name, 'N(\d+)R', 'tokens');
    
    if isempty(tokens)
        warning('Cannot extract N from %s, skipping.', files(k).name);
        continue
    end
    
    Nsamples = str2double(tokens{1}{1});

    % Add NumSamples column
    Ncol = repmat(Nsamples, height(T), 1);
    T = addvars(T, Ncol, 'After', 'Method', ...
        'NewVariableNames', col_name);
    
    % Append
    Tall = [Tall; T];
end

% Sort table
Tall = sortrows(Tall, col_name);

%% =========================================================
%        PROCESS LQR ROWS (MEAN AGGREGATION)
%% =========================================================

% Identify LQR rows
idx = strcmp(Tall.Method, 'Nominal point LQR');
LQR = Tall(idx, :);

% ----------- Helper functions -----------
% Extract numeric from string
toNum = @(c) str2double(c);

% Extract mean and std from "a ± b"
parse_pm = @(cellArray) cellfun(@(s) ...
    regexp(s, '([0-9.eE+-]+)\s*±\s*([0-9.eE+-]+)', 'tokens', 'once'), ...
    cellArray, 'UniformOutput', false);

% ----------- Simple numeric columns -----------
MSE_nom = toNum(LQR.MSE_traj_nominal_rho);
MSE_rand = toNum(LQR.MSE_traj_random_rho);
time = toNum(LQR.("Solver time [sec]"));

% ----------- Parse ± columns -----------
tokens_nom = parse_pm(LQR.SS_error_nominal_rho);
tokens_rand = parse_pm(LQR.SS_error_random_rho);

means_nom = cellfun(@(t) str2double(t{1}), tokens_nom);
stds_nom  = cellfun(@(t) str2double(t{2}), tokens_nom);

means_rand = cellfun(@(t) str2double(t{1}), tokens_rand);
stds_rand  = cellfun(@(t) str2double(t{2}), tokens_rand);

% ----------- Aggregate -----------
mean_table = table( ...
    {'Nominal point LQR (mean)'}, ...
    0, ...
    {sprintf('%.2e', mean(MSE_nom))}, ...
    {sprintf('%.2e', mean(MSE_rand))}, ...
    {sprintf('%.2e ± %.2e', mean(means_nom), mean(stds_nom))}, ...
    {sprintf('%.2e ± %.2e', mean(means_rand), mean(stds_rand))}, ...
    {sprintf('%.3f', mean(time))}, ...
    'VariableNames', Tall.Properties.VariableNames);

%% =========================================================
%        FINAL TABLE: REMOVE OLD LQR + ADD MEAN ROW
%% =========================================================

Tall_clean = Tall(~idx, :);         % remove all LQR rows
Tall_final = [Tall_clean; mean_table];  % append mean LQR at end

% Display final table
disp(Tall_final)
Tall = Tall_final;
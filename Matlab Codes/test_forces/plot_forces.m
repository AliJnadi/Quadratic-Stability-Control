%% ================= MERGE & PLOT RESULTS ===================
clc; clear; close all;

Font_Size = 14;

n_states = 16;     % length of force selector mask
n_files  = 15;     % F1 ... F15

results = [];      % master struct
idx = 1;

%% ================= LOAD ALL MAT FILES =====================
for i = 1:n_files
    fname = sprintf('F%d.mat', i);
    if ~isfile(fname)
        warning('%s not found, skipping.', fname);
        continue;
    end
    
    % Extract table T
    S = load(fname);
    if ~isfield(S, 'T')
        warning('No table T in %s, skipping.', files(k).name);
        continue
    end
    T = S.T;
    
    % --- Extract force selector
    bin_r = dec2bin(i, 4);
    fs = double(bin_r) - '0';

    % --- Count number of actuators
    numAct = sum(fs == 1);

    % ---------- QSC ----------
    results(idx).method = 'QSC';
    results(idx).numAct = numAct;
    results(idx).MSE    = str2double(T{1,3});
    results(idx).SS     = parse_mean(T{1,5});
    idx = idx + 1;

    % ---------- LQR ----------
    results(idx).method = 'LQR';
    results(idx).numAct = numAct;
    results(idx).MSE    = str2double(T{2,3});
    results(idx).SS     = parse_mean(T{2,5});
    idx = idx + 1;
end

fprintf('Loaded %d method entries\n', numel(results));

%% ================= AGGREGATE BY ACTUATOR COUNT =============
actCounts = 1:4;

MSE_QSC = nan(1,4);
MSE_LQR = nan(1,4);
SS_QSC  = nan(1,4);
SS_LQR  = nan(1,4);

for k = actCounts
    idxQ = ([results.numAct] == k) & strcmp({results.method},'QSC');
    idxL = ([results.numAct] == k) & strcmp({results.method},'LQR');

    if any(idxQ)
        MSE_QSC(k) = mean([results(idxQ).MSE]);
        SS_QSC(k)  = mean([results(idxQ).SS]);
    end
    if any(idxL)
        MSE_LQR(k) = mean([results(idxL).MSE]);
        SS_LQR(k)  = mean([results(idxL).SS]);
    end
end

%% ===================== PLOT: MSE ==========================
figure('Color','w','Position',[100 100 600 420]);
hold on; grid on;

ax1 = gca();
ax1.FontName = 'Times New Roman';
ax1.FontSize = Font_Size;

bar(actCounts, [MSE_LQR' MSE_QSC'], 'grouped');
xlabel('Number of actuators','FontSize',Font_Size);
ylabel('Trajectory MSE [m^2] (random \rho)','FontSize',12);

legend('LQR','QSC', ...
       'Location','northoutside', ...
       'NumColumns', 2);

ax1.Legend.Title.String = 'Trajectory MSE vs Actuator Count';
ax1.Legend.Title.FontSize = Font_Size;

print('-dpng', '-painters', 'trajectory_MSE_vs_actuator_count.png');
print('-depsc2', '-painters', 'trajectory_MSE_vs_actuator_count.eps');

% title('Trajectory MSE vs Actuator Count','FontSize',14);
%% ================== PLOT: SS ERROR ========================
figure('Color','w','Position',[100 100 600 420]);
hold on; grid on;

ax1 = gca();
ax1.FontName = 'Times New Roman';
ax1.FontSize = Font_Size;

bar(actCounts, [SS_LQR' SS_QSC'], 'grouped');
xlabel('Number of actuators', ...
       'FontName', 'Times New Roman', ...
       'FontSize', Font_Size);
ylabel('Steady-state error [m] (random \rho)', ...
       'FontName', 'Times New Roman', ...
       'FontSize', Font_Size)

ax1.LineWidth = 1;

legend('LQR','QSC', ...
       'Location','northoutside', ...
       'NumColumns', 2);
   
ax1.Legend.Title.String = 'Steady-State Error vs Actuator Count';
ax1.Legend.Title.FontSize = Font_Size;


print('-dpng', '-painters', 'steady_state_error_vs_actuator_count.png');
print('-depsc2', '-painters', 'steady_state_error_vs_actuator_count.eps');

% title('Steady-State Error vs Actuator Count','FontSize',14);

%% ===================== DONE ===============================
disp('Plots generated successfully.');

Tall = table(MSE_QSC', MSE_LQR', SS_QSC', SS_LQR', 'VariableNames', ["MSE_QSC", "MSE_LQR", "SS_QSC", "SS_LQR"]);
Tall = addvars(Tall, (1:height(Tall))', 'Before', 'MSE_QSC', 'NewVariableNames', 'RowNumber');
writetable(Tall, 'forces.csv');
%% ============== HELPER FUNCTION ===========================
function m = parse_mean(cellEntry)
% Extract mean from string like '1.37e-01 ± 8.25e-02'
    if isnumeric(cellEntry)
        m = cellEntry;
        return;
    end
    s = char(cellEntry);
    tokens = regexp(s,'^([0-9eE\.\+\-]+)','tokens');
    m = str2double(tokens{1});
end

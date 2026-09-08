Font_Size = 14;
%% ============================================================
%  Load and prepare data
% ============================================================
merge_results;   % produces Tall

%% -----------------------------
% Convert string columns to numeric
% -----------------------------
Tall.MSE_nom_val = cellfun(@str2double, Tall.MSE_traj_nominal_rho);
Tall.MSE_ran_val = cellfun(@str2double, Tall.MSE_traj_random_rho);

extractMean = @(s) sscanf(s, '%e ± %*e');

Tall.SS_nom_mean = cellfun(extractMean, Tall.SS_error_nominal_rho);
Tall.SS_ran_mean = cellfun(extractMean, Tall.SS_error_random_rho);

Tall.SolverTime = cellfun(@str2double, Tall.('Solver time [sec]'));

methods = unique(Tall.Method);
%% -----------------------------
% Separate QSC and LQR
% -----------------------------
isLQR = contains(Tall.Method, 'LQR');

Tall_QSC = Tall(~isLQR, :);
Tall_LQR = Tall(isLQR, :);

Tall_QSC = Tall_QSC(Tall_QSC.NumSamples > 0, :);

Ns = Tall_QSC.NumSamples;
xrange = [min(Ns) max(Ns)];

idx = strcmp(Tall_QSC.Method, 'Quadratic stability control');

%% ============================================================
% (1) MSE vs Number of Samples
% ============================================================
figure('Color','w','Position',[100 100 600 420]);
hold on; grid on;

h = [];
vals = [];

% QSC nominal
h(end+1) = plot(Ns(idx), Tall_QSC.MSE_nom_val(idx), '-o', 'LineWidth', 2);
vals(end+1) = max(Tall_QSC.MSE_nom_val(idx));

% QSC random
h(end+1) = plot(Ns(idx), Tall_QSC.MSE_ran_val(idx), '--s', 'LineWidth', 2);
vals(end+1) = max(Tall_QSC.MSE_ran_val(idx));

% LQR nominal
y_nom = Tall_LQR.MSE_nom_val(1);
h(end+1) = plot(xrange, [y_nom y_nom], '-', 'LineWidth', 2);
vals(end+1) = y_nom;

% LQR random
y_ran = Tall_LQR.MSE_ran_val(1);
h(end+1) = plot(xrange, [y_ran y_ran], '--', 'LineWidth', 2);
vals(end+1) = y_ran;

% Sort legend
[~, order] = sort(vals, 'descend');

labels = { ...
    'QSC (proposed, nominal)', ...
    'QSC (proposed, random)', ...
    'LQR (baseline, nominal)', ...
    'LQR (baseline, random)'};

set(gca,'YScale','log')

ax1 = gca();
ax1.FontName = 'Times New Roman';
ax1.FontSize = Font_Size;

legend(h(order), labels(order), 'Location','best');
   
ax1.Legend.Title.String = 'Trajectory MSE vs Number of Samples';
ax1.Legend.Title.FontSize = Font_Size;

xlabel('Number of samples N','FontName','Times New Roman','FontSize',Font_Size);
ylabel('Trajectory MSE [m^2] (log scale)','FontName','Times New Roman','FontSize',Font_Size);

ax1.LineWidth = 1;
grid on;

set(gca,'FontName','Times New Roman','FontSize',Font_Size,'LineWidth',1);

savefig('trajectory_mse_vs_samples.fig');
print('-dpng', '-painters', 'trajectory_mse_vs_samples.png');
print('-depsc2', '-painters', 'trajectory_mse_vs_samples.eps');

%% ============================================================
% (2) Steady-State Error
% ============================================================
figure('Color','w','Position',[100 100 600 420]);
hold on; grid on;

h = [];
vals = [];

% QSC
h(end+1) = plot(Ns(idx), Tall_QSC.SS_ran_mean(idx), '-o', 'LineWidth', 2);
vals(end+1) = max(Tall_QSC.SS_ran_mean(idx));

% LQR
y_lqr = Tall_LQR.SS_ran_mean(1);
h(end+1) = plot(xrange, [y_lqr y_lqr], '-', 'LineWidth', 2);
vals(end+1) = y_lqr;

[~, order] = sort(vals, 'descend');

labels = {'QSC (proposed)', 'LQR (baseline)'};
set(gca,'YScale','log')

ax1 = gca();
ax1.FontName = 'Times New Roman';
ax1.FontSize = Font_Size;

legend(h(order), labels(order), 'Location','best');

x1.Legend.Title.String = 'SSE under Random Uncertainty';
ax1.Legend.Title.FontSize = Font_Size;
 

xlabel('Number of samples N','FontName','Times New Roman','FontSize',Font_Size);
ylabel('SSE [m] (log scale)','FontName','Times New Roman','FontSize',Font_Size);

set(gca,'FontName','Times New Roman','FontSize',Font_Size,'LineWidth',1);

ax1.LineWidth = 1;
grid on;

savefig('ss_error_vs_samples.fig');
print('-dpng', '-painters', 'ss_error_vs_samples.png');
print('-depsc2', '-painters', 'ss_error_vs_samples.eps');

%% ============================================================
% (3) Solver Time
% ============================================================
figure('Color','w','Position',[100 100 600 420]);
hold on; grid on;

h = [];
vals = [];

% QSC
h(end+1) = plot(Ns(idx), Tall_QSC.SolverTime(idx), '-o', 'LineWidth', 2);
vals(end+1) = max(Tall_QSC.SolverTime(idx));

% LQR
y_time = Tall_LQR.SolverTime(1);
h(end+1) = plot(xrange, [y_time y_time], '-', 'LineWidth', 2);
vals(end+1) = y_time;

[~, order] = sort(vals, 'descend');

labels = {'QSC (proposed)', 'LQR (baseline)'};

ax1 = gca();
ax1.FontName = 'Times New Roman';
ax1.FontSize = Font_Size;

legend(h(order), labels(order), 'Location','best');

ax1.Legend.Title.String = 'Solver Time vs Number of Samples';
ax1.Legend.Title.FontSize = Font_Size;

ax1.LineWidth = 1;
grid on;

xlabel('Number of samples N','FontName','Times New Roman','FontSize',Font_Size);
ylabel('Solver time [sec]','FontName','Times New Roman','FontSize',Font_Size);

set(gca,'FontName','Times New Roman','FontSize',Font_Size,'LineWidth',1);

savefig('time_vs_samples.fig');
print('-dpng', '-painters', 'time_vs_samples.png');
print('-depsc2', '-painters', 'time_vs_samples.eps');
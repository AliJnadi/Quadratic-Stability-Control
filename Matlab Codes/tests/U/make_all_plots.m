Font_Size = 14;

%% ============================================================
%  Separate figures using PRINT
% ============================================================
merge_results_U;

%% -----------------------------
% Convert string columns to numeric
% -----------------------------

% MSE
Tall.MSE_nom_val = cellfun(@str2double, Tall.MSE_traj_nominal_rho);
Tall.MSE_ran_val = cellfun(@str2double, Tall.MSE_traj_random_rho);

% Steady-state mean extraction: "x.xx ± y.yy"
extractMean = @(s) sscanf(s, '%e ± %*e');
Tall.SS_nom_mean = cellfun(extractMean, Tall.SS_error_nominal_rho);
Tall.SS_ran_mean = cellfun(extractMean, Tall.SS_error_random_rho);

% Solver time
Tall.SolverTime = cellfun(@str2double, Tall.('Solver time [sec]'));

methods = unique(Tall.Method);

%% -----------------------------
% Separate QSC and LQR
% -----------------------------
isLQR = contains(Tall.Method, 'LQR');

Tall_QSC = Tall(~isLQR, :);
Tall_LQR = Tall(isLQR, :);

R = Tall_QSC.SamplingRadius;
xrange = [min(R) max(R)];

LQR_idx = strcmp(Tall_LQR.Method, 'Nominal point LQR');
QSC_idx = strcmp(Tall_QSC.Method, 'Quadratic stability control');

%% ============================================================
% (1) MSE vs Number of Samples
% ============================================================
figure('Color','w','Position',[100 100 600 420]);
hold on; grid on;

h = [];
vals = [];

% QSC nominal
y_nom =  mean(Tall_QSC.MSE_nom_val(QSC_idx));
h(end+1) = plot(xrange, [y_nom y_nom], '-', 'LineWidth', 2);
vals(end+1) = y_nom;

% QSC random
h(end+1) = plot(R(QSC_idx), Tall_QSC.MSE_ran_val(QSC_idx), '-s', 'LineWidth', 2);
vals(end+1) = max(Tall_QSC.MSE_ran_val(QSC_idx));


% LQR nominal
y_nom =  mean(Tall_LQR.MSE_nom_val(LQR_idx));
h(end+1) = plot(xrange, [y_nom y_nom], '--', 'LineWidth', 2);
vals(end+1) = y_nom;

% LQR random
h(end+1) = plot(R(LQR_idx), Tall_LQR.MSE_ran_val(LQR_idx), '--s', 'LineWidth', 2);
vals(end+1) = max(Tall_LQR.MSE_ran_val(LQR_idx));

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
   

ax1.Legend.Title.String = 'Trajectory MSE vs Parameter Uncertainity';
ax1.Legend.Title.FontSize = Font_Size;

xlabel('Parameter Uncertainity U [%]', ...
       'FontName', 'Times New Roman', ...
       'FontSize', Font_Size);
   
ylabel('Trajectory MSE [m^2] (log scale)', ...
       'FontName', 'Times New Roman', ...
       'FontSize', Font_Size);  

xlim(xrange);
ax1.LineWidth = 1;
grid on;

set(gca,'FontName','Times New Roman','FontSize',Font_Size,'LineWidth',1);

savefig('trajectory_mse_vs_uncertinity.fig');
print('-dpng', '-painters', 'trajectory_mse_vs_uncertinity.png');
print('-depsc2', '-painters', 'trajectory_mse_vs_uncertinity.eps');
%% ============================================================
% (2) Steady-State Error vs Number of Samples
% ============================================================
figure('Color','w','Position',[100 100 600 420]);
hold on; grid on;

h = [];
vals = [];

% QSC
h(end+1) = plot(R(QSC_idx), Tall_QSC.SS_ran_mean(QSC_idx), '-o', 'LineWidth', 2);
vals(end+1) = max(Tall_QSC.SS_ran_mean(QSC_idx));

% LQR
h(end+1) = plot(R(LQR_idx), Tall_LQR.SS_ran_mean(LQR_idx), '-o', 'LineWidth', 2);
vals(end+1) = max(Tall_LQR.SS_ran_mean(LQR_idx));

ax1 = gca();
ax1.FontName = 'Times New Roman';
ax1.FontSize = Font_Size;

[~, order] = sort(vals, 'descend');

labels = {'QSC (proposed)', 'LQR (baseline)'};
legend(h(order), labels(order), 'Location','best');

set(gca,'YScale','log')

ax1.Legend.Title.String = 'SSE under Random Uncertainty';
ax1.Legend.Title.FontSize = Font_Size;

xlabel('Parameter Uncertainity U [%]', ...
       'FontName', 'Times New Roman', ...
       'FontSize', Font_Size);
   
ylabel('SSE [m] (log scale)', ...
       'FontName', 'Times New Roman', ...
       'FontSize', Font_Size);

xlim(xrange);
ax1.LineWidth = 1;
grid on;

set(gca,'FontName','Times New Roman','FontSize',Font_Size,'LineWidth',1);

savefig('ss_error_vs_uncertinity.fig');
print('-dpng', '-painters', 'ss_error_vs_uncertinity.png');
print('-depsc2', '-painters', 'ss_error_vs_uncertinity.eps');
%% ============================================================
% (3) Solver Time vs Number of Samples
% ============================================================
figure('Color','w','Position',[100 100 600 420]);
hold on;

h = [];
vals = [];

% QSC
y_time = mean(Tall_QSC.SolverTime(QSC_idx));
h(end+1) = plot(xrange, [y_time y_time], '-', 'LineWidth', 2);
vals(end+1) = y_time

% LQR
y_time = mean(Tall_LQR.SolverTime(LQR_idx));
h(end+1) = plot(xrange, [y_time y_time], '-', 'LineWidth', 2);
vals(end+1) = y_time

[~, order] = sort(vals, 'descend');

labels = {'QSC (proposed)', 'LQR (baseline)'};

ax1 = gca();
ax1.FontName = 'Times New Roman';
ax1.FontSize = Font_Size;

legend(h(order), labels(order), 'Location','best');

ax1.Legend.Title.String = 'Solver Time vs Parameter Uncertainity';
ax1.Legend.Title.FontSize = Font_Size;

xlabel('Parameter Uncertainity U [%]', ...
       'FontName', 'Times New Roman', ...
       'FontSize', Font_Size);
   
ylabel('Solver time [sec]', ...
       'FontName', 'Times New Roman', ...
       'FontSize', Font_Size);

xlim(xrange);
ax1.LineWidth = 1;
grid on;

set(gca,'FontName','Times New Roman','FontSize',Font_Size,'LineWidth',1);

savefig('time_vs_uncertinity.fig');
print('-dpng', '-painters', 'time_vs_uncertinity.png');
print('-depsc2', '-painters', 'time_vs_uncertinity.eps');
% title('Solver Time vs Number of Samples','Interpreter','latex');

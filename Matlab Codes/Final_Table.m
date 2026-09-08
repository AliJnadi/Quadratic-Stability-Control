Method = {
    'Quadratic stability control'
    'Nominal point LQR'
};

MSE_nom = {
    sprintf('%.2e', stats_QSC_nom.traj_mse_mean)
    sprintf('%.2e', stats_LQR_nom.traj_mse_mean)
};

MSE_ran = {
    sprintf('%.2e', stats_QSC_ran.traj_mse_mean)
    sprintf('%.2e', stats_LQR_ran.traj_mse_mean)
};

SS_nom = {
    sprintf('%.2e ± %.2e', stats_QSC_nom.ss_mean,  stats_QSC_nom.ss_std)
    sprintf('%.2e ± %.2e', stats_LQR_nom.ss_mean,  stats_LQR_nom.ss_std)
};

SS_ran = {
    sprintf('%.2e ± %.2e', stats_QSC_ran.ss_mean,  stats_QSC_ran.ss_std)
    sprintf('%.2e ± %.2e', stats_LQR_ran.ss_mean,  stats_LQR_ran.ss_std)  
};

Solver_time = {
    sprintf('%.3f', info_QSC.elapsed_time)
    sprintf('%.3f', info_LQR.elapsed_time)
};

T = table( ...
    Method, ...
    MSE_nom, ...
    MSE_ran, ...
    SS_nom, ...
    SS_ran, ...
    Solver_time, ...
    'VariableNames', { ...
        'Method', ...
        'MSE_traj_nominal_rho', ...
        'MSE_traj_random_rho', ...
        'SS_error_nominal_rho', ...
        'SS_error_random_rho', ...
        'Solver time [sec]' ...
    });

display(T)
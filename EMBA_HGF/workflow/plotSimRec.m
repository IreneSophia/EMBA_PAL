function [] = plotSimRec(modSpace, nSim, saveDir)
% 'sim_rec_avg_logRT', 'sim_rec_logRT_20sub_model': plot RT trajectories for 
% all simulatd subjecta as well as the average
% Based on hessetal_spirl_analysis/job_runner_plot_results.m, lines
% 855-917%
% 
% !! CAUTION !!: This is HARDCODED for the models used in the EMBA project. 
%
% INPUT
%   modSpace            struct              model space
%
%   nSim                integer             number of simulated subjects
%
%   saveDir             char array          base output directory
%-----------------------------------------------------------------------------
%
% Copyright (C) 2024 Anna Yurova, Irene Sophia Plank, LMU University Hospital;
%               based on the hessetal_spirl_analysis toolbox by Alex Hess (2024), TNU, ETHZ
%
%
% This file is part of the EMBA_HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. For further details, see the file LICENSE or <http://www.gnu.org/licenses/>.
%-----------------------------------------------------------------------------

nModels = size(modSpace, 2);

res = struct();
res.rec = load(fullfile(saveDir, 'sim', ['recovery_analysis']));

nTrials = numel(res.rec.est(1,1,1).data.y);
trials = 1:nTrials;

input = res.rec.est(1,1,1).data.u(:,1);
input(input == 1) = 7.5;
input(input == 0) = 5.5;


%% REC: plot AVG rt trajectories & fits (same model sim + est)
figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
for m = 1:nModels
    for n = 1:nSim
        y_sim_mat(:,n) = res.rec.est(m,n,m).data.y;
        yhat_mat(:,n) = res.rec.est(m,n,m).data.optim.yhat;
    end
    
    mean_logRT_sim = mean(y_sim_mat, 2, 'omitnan');
    var_logRT_sim = var(y_sim_mat, 0, 2, 'omitnan');
    logupper_sim = mean_logRT_sim + sqrt(var_logRT_sim);
    loglower_sim = mean_logRT_sim - sqrt(var_logRT_sim);
    mean_yhat = mean(yhat_mat, 2, 'omitnan');
    var_yhat = var(yhat_mat, 0, 2, 'omitnan');
    upper_yhat = mean_yhat + sqrt(var_yhat);
    lower_yhat = mean_yhat - sqrt(var_yhat);
    
    subplot(2,2,m)
    plot(trials, mean_logRT_sim, 'LineWidth', 2, 'color', 'b'); %, 'color', options.col.tnub)
    hold on
    fill([trials, fliplr(trials)], [(logupper_sim)', fliplr((loglower_sim)')], ...
             'b', 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
    plot(trials, mean_yhat, 'LineWidth', 2)
    fill([trials, fliplr(trials)], [(upper_yhat)', fliplr((lower_yhat)')], ...
             'r', 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
    % add the input and labels and such
    scatter(1:length(input), input, 5, 'o','MarkerEdgeColor','k','MarkerEdgeAlpha', 0.6, 'LineWidth', 0.1)
    xlim([1 nTrials])
    ylim([5.25 7.75])
    ylabel('logRT [ms]', 'FontSize', 14, 'Interpreter','latex')
    xlabel('trials', 'FontSize', 14, 'Interpreter','latex')
    txt = "model " + modSpace(m).name;
    title(txt, 'FontSize', 20, 'Interpreter','latex')
end

legend('$log(y_{rt})$', '$sd(log(y_{rt}))$', '$log(\hat{y}_{rt})$', '$sd(log(\hat{y}_{rt}))$', 'Interpreter','latex', 'Position', [0.94 0.48 0.03 0.07], 'FontSize', 14)
figdir = fullfile(saveDir, 'figures', 'sim_rec');
if ~exist(figdir, 'dir')
    mkdir(figdir)
end
print(strcat(figdir, filesep, 'sim_rec_avg_logRT'), '-dpng');
print(strcat(figdir, filesep, 'sim_rec_avg_logRT'), '-dsvg');
close;


end

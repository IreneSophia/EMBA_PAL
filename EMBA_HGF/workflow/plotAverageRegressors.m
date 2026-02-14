function [] = plotAverageRegressors(modSpace, subDir, data, saveDir)
% Plot average regressors: 'avg_logRT_regressors_mod', 
% Based on hessetal_spirl_analysis/job_runner_plot_results.m, lines
% 178-191,215-248
%
% !! CAUTION !!: This is HARDCODED for the models used in the EMBA project. 
% 
% Main changes:
%    - Changed the formula for all outputted values of the 'regressors model' plot 
%     to the ones that I think correspond to our model. In particular, I
%     took the formulas from the tapas_logrt_linear_binary.m file
%
%
% INPUT
%   modSpace            struct              model space
%
%   subDir              char array          subdirectory containing results
%
%   data                struct              data (m x nSubs)
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

load(fullfile(saveDir, subDir, 'full_results'), 'res');

nSubs   = size(res.est, 2);
nTrials = size(data(1).u,1);
trials = 1:nTrials;
logRT_mat = zeros(nTrials, nSubs);

% get the input sequence
u = data(end).u(:,1);

for i = 1:nTrials
    for n = 1:nSubs
        logRT_mat(i, n) = data(n).y(i);
    end
end

mean_logRT = mean(logRT_mat, 2, 'omitnan');
var_logRT  = var(logRT_mat, 0, 2, 'omitnan');
logupper = mean_logRT + sqrt(var_logRT);
loglower = mean_logRT - sqrt(var_logRT);

for m = 1:nModels
    yhat_mat = nan(nTrials, nSubs);
    for n = 1:nSubs
        yhat_mat(:,n) = res.est(m,n).optim.yhat;
    end
    mod(m).mean_yhat  = mean(yhat_mat, 2, 'omitnan');
    mod(m).var_yhat   = var(yhat_mat, 0, 2, 'omitnan');
    mod(m).upper_yhat = mod(m).mean_yhat + sqrt(mod(m).var_yhat);
    mod(m).lower_yhat = mod(m).mean_yhat - sqrt(mod(m).var_yhat);

end

for m = 1:nModels
    %% 'avg log rt fits' plot
    figure
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.3 0 0.45 1]);
    subplot(2,1,1)
    plot(1:nTrials, mean_logRT, 'LineWidth', 2); %, 'color', options.col.tnub)
    hold on
    fill([trials, fliplr(trials)], [(logupper)', fliplr((loglower)')], ...
             'b', 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
    plot(trials, mod(m).mean_yhat, 'LineWidth', 2)
    fill([trials, fliplr(trials)], [(mod(m).upper_yhat)', fliplr((mod(m).lower_yhat)')], ...
             'r', 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
    % input
    input = u(:,1)*6.6;
    input(input == 0) = 6.0;
    scatter(1:length(input), input, 5, 'o','MarkerEdgeColor','k','MarkerEdgeAlpha', 0.6, 'LineWidth', 0.1)
    % add legends and stuff
    legend('$log(y_{rt})$', '$sd(log(y_{rt}))$', '$log(\hat{y}_{rt})$', '$sd(log(\hat{y}_{rt}))$', ...
        'Interpreter','latex', 'Location', 'southoutside', 'Orientation','horizontal')
    xlim([1 nTrials])
    ylim([5.8 6.9])
    ylabel('logRT [ms]', 'Interpreter','latex')
    title('avg log rt fits', 'Interpreter','latex', 'Fontsize', 20)

    %% 'regressors model ' plot: HARDCODED FOR PARTICULAR MODELS
    for n = 1:nSubs

        % formulas from tapas_logrt_linear_binary.m
        m1hreg = res.est(m,n).traj.muhat(:,1);
        poo    = m1hreg.^res.est(m,n).u(:,1).*(1-m1hreg).^(1-res.est(m,n).u(:,1)); % probability of observed outcome
        surp_mat(:,n)   = -log2(poo);

        unc1_mat(:,n) = res.est(m,n).traj.sahat(:,1);

        mu2 = res.est(m,n).traj.mu(:,2);
        sa2 = res.est(m,n).traj.sa(:,2);
        unc2_mat(:,n) = tapas_sgm(mu2, 1).*(1 -tapas_sgm(mu2, 1)).*sa2;

        if modSpace(m).obs == "emba_logrt_linear_binary_C"
            pwpe_mat(:,n) = res.est(m,n).traj.epsi(:,3);
        end

        if modSpace(m).obs == "tapas_logrt_linear_binary"
            mu3 = res.est(m,n).traj.mu(:,3);
            vol_mat(:,n)  = tapas_sgm(mu2, 1).*(1-tapas_sgm(mu2, 1)).*exp(mu3);
        end
        
    end

    % AVG regressors
    subplot(2,1,2)
    hold on

    % Changed the legend too, to correspond to the Lawson model
    plot(trials, rescale(mean(surp_mat, 2, 'omitnan'),0,1), 'Color', '#648FFF');    
    plot(trials, rescale(mean(unc1_mat, 2, 'omitnan'),0,1), 'Color', '#785EF0');
    plot(trials, rescale(mean(unc2_mat, 2, 'omitnan'),0,1), 'Color', '#DC267F');
    
    % legend text
    txt = {'$surprise$', '$unc_1$', '$unc_2$'};

    if modSpace(m).obs == "emba_logrt_linear_binary_C"
        % plot precision-weighted prediction error
        plot(trials, rescale(mean(pwpe_mat, 2, 'omitnan'),0,1), 'Color', '#FE6100');
        % add to legend text
        txt = [txt, '$pwPE$'];
    end


    if modSpace(m).obs == "tapas_logrt_linear_binary"
        % plot phasic volatility
        plot(trials, rescale(mean(vol_mat, 2, 'omitnan'),0,1), 'Color', '#FFB000')
        % add to legend text
        txt = [txt, '$volatility$'];
    end
    % input
    input = u(:,1)*1.2;
    input(u == 0) = -0.2;
    scatter(1:length(input), input, 5, 'o','MarkerEdgeColor','k','MarkerEdgeAlpha', 0.6, 'LineWidth', 0.1)
    % add legend
    legend(txt, 'Interpreter', 'Latex', 'Location', 'southoutside', 'Orientation','horizontal');
    % specify limits    
    ylim([-0.3, 1.3])
    % add labels and stuff
    ylabel('arbitrary units', 'Interpreter','latex')
    xlabel('trials', 'Interpreter','latex')
    txt = 'response model regressors';
    title(txt, 'Interpreter','latex', 'Fontsize', 20)
    
    figdir = fullfile(saveDir, 'figures', subDir, 'regressors');
    if ~exist(figdir, 'dir')
        mkdir(figdir)
    end
    print(strcat(figdir,filesep,'avg_logRT_regressors_mod_',modSpace(m).name), '-dpng');
    print(strcat(figdir,filesep,'avg_logRT_regressors_mod_',modSpace(m).name), '-dsvg');
    close;
end

end

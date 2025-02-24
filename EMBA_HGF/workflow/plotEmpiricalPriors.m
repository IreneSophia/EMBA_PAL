function [] = plotEmpiricalPriors(modSpace, saveDir)
% Plot empirical priors
% Based on hessetal_spirl_analysis/job_runner_paper_figs.m, lines
% 404 - 579 -> produces figure 5A in their paper
%
% !! CAUTION !!: This is HARDCODED for the models used in the EMBA project. 
%
% INPUT
%   modSpace            struct              model space
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
empiricalPriors = load(fullfile(saveDir, 'pilots', 'pilot_priors'));
model = empiricalPriors.mod;

%% FIG: empirical prior distributions

figpath = fullfile(saveDir, 'figures', 'priors');
if ~exist(figpath, 'dir')
    mkdir(figpath);
end

% x grid
x_min = -50;
x_max = 50;
x = x_min:0.001:x_max;

for m = 1:size(modSpace, 2)

    % get pdf of prc model params
    for j = 1:size(modSpace(m).prc_idx,2)
        prc(j).y = normpdf(x, model(m).prc_robmean(1,j),...
            sqrt(model(m).prc_robvar(1,j)));
        prc(j).y_prior = normpdf(x, model(m).prc_priormus(j),...
            sqrt(model(m).prc_priorsas(j)));
    end
    
    % get pdf of obs model params
    for k = 1:size(modSpace(m).obs_idx,2)
        obs(k).y = normpdf(x, model(m).obs_robmean(1,k),...
            sqrt(model(m).obs_robvar(1,k)));
        obs(k).y_prior = normpdf(x, model(m).obs_priormus(k),...
            sqrt(model(m).obs_priorsas(k)));
    end
    
    figure
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    
    %% PERCEPTION MODEL
    
    % plot om2
    subplot(3,3,1)
    j = 1;
    plot(x, prc(j).y, 'k', 'LineWidth', 1)
    hold on
    plot(x, prc(j).y_prior, 'k--')
    % add pilot participants
    plot(model(m).prc_est(:,j), -0.05, 'k.')
    hold off
    %ylim([-0.1 0.25])
    xlim([-15 5])
    t = '$\omega_2$';
    title(t, 'interpreter','latex', 'FontSize', 20);

    % add a legend
    legend('empirical priors', 'default priors',...
        'Interpreter','latex', 'Position', [0.05 0.005 0.9 0.05], 'Orientation', 'horizontal', 'FontSize', 14)
    
    % plot om3
    if modSpace(m).prc_config.n_levels > 2
        subplot(3,3,2);
        j = 2;
        plot(x, prc(j).y, 'k', 'LineWidth', 1)
        hold on
        plot(x, prc(j).y_prior, 'k--')
        plot(model(m).prc_est(:,j), -0.05, 'k.')
        hold off
        %ylim([-0.1 0.25])
        xlim([-15 5])
        t = '$\omega_3$';
        title(t, 'interpreter','latex', 'FontSize', 20);
    end

    % plot alpha
    if length(modSpace(m).prc_idx) > 2
        subplot(3,3,3);
        j = 3;
        plot(x, prc(j).y, 'k', 'LineWidth', 1)
        hold on
        plot(x, prc(j).y_prior, 'k--')
        plot(model(m).prc_est(:,j), -0.05, 'k.')
        hold off
        %ylim([-0.1 0.25])
        xlim([-4 2])
        t = '$log(\alpha)$';
        title(t, 'interpreter','latex', 'FontSize', 20);

    end

    %% RESPONSE MODEL
    
    % plot beta0 (intercept)
    subplot(3,3,4)
    k = 1;
    plot(x, obs(k).y, 'k', 'LineWidth', 1)
    hold on
    plot(x, obs(k).y_prior, 'k--')
    plot(model(m).obs_est(:,k), -0.05, 'k.')
    hold off
    %ylim([-0.1 2])
    xlim([log(500)-2 log(500)+2])
    t = '$\beta_0$';
    title(t, 'interpreter','latex', 'FontSize', 20);
    
    % plot beta1 (Shannon surprise)
    subplot(3,3,5)
    k = 2;
    plot(x, obs(k).y, 'k', 'LineWidth', 1)
    hold on
    plot(x, obs(k).y_prior, 'k--')
    plot(model(m).obs_est(:,k), -0.05, 'k.')
    hold off
    %ylim([-0.1 16])
    xlim([-0.25 0.25])
    t = '$\beta_{surprise}$';
    title(t, 'interpreter','latex', 'FontSize', 20);

    if length(modSpace(m).obs_config.priormus) == 6
    
        % plot beta2 
        subplot(3,3,6)
        k = 3;
        plot(x, obs(k).y, 'k', 'LineWidth', 1)
        hold on
        plot(x, obs(k).y_prior, 'k--')
        plot(model(m).obs_est(:,k), -0.05, 'k.')
        hold off
        %ylim([-0.1 0.4])
        xlim([-6 6])
        t = '$\beta_{unc1}$';
        title(t, 'interpreter','latex', 'FontSize', 20);
        
        % plot beta3 
        subplot(3,3,7)
        k = 4;
        plot(x, obs(k).y, 'k', 'LineWidth', 1)
        hold on
        plot(x, obs(k).y_prior, 'k--')
        plot(model(m).obs_est(:,k), -0.05, 'k.')
        hold off
        %ylim([-0.1 0.3])
        xlim([-8 8])
        t = '$\beta_{unc2}$';
        title(t, 'interpreter','latex', 'FontSize', 20);
        k = 5;

    else

        % plot beta pwPE 
        subplot(3,3,6)
        k = 3;
        plot(x, obs(k).y, 'k', 'LineWidth', 1)
        hold on
        plot(x, obs(k).y_prior, 'k--')
        plot(model(m).obs_est(:,k), -0.05, 'k.')
        hold off
        %ylim([-0.1 0.4])
        xlim([-6 6])
        t = '$\beta_{pwPE}$';
        title(t, 'interpreter','latex', 'FontSize', 20);
        k = 4;

    end

    % plot phasic volatility 
    subplot(3,3,8)
    plot(x, obs(k).y, 'k', 'LineWidth', 1)
    hold on
    plot(x, obs(k).y_prior, 'k--')
    plot(model(m).obs_est(:,k), -0.05, 'k.')
    hold off
    %ylim([-0.1 0.7])
    xlim([-4 4])
    t = '$\beta_{volatility}$';
    title(t, 'interpreter','latex', 'FontSize', 20);
    
    % plot Sigma
    subplot(3,3,9)
    % Set index depending on the number of parameters in the corresponding
    % model
    k = length(model(m).obs_robmean);
    plot(x, obs(k).y, 'k', 'LineWidth', 1)
    hold on
    plot(x, obs(k).y_prior, 'k--')
    plot(model(m).obs_est(:,k), -0.05, 'k.')
    hold off
    %ylim([-0.1 2])
    xlim([-5 4])
    t = '$log(\Sigma)$';
    title(t, 'interpreter','latex', 'FontSize', 20);
    
    print([figpath filesep char(strcat('fig_empirical_priors_', modSpace(m).name))], '-dpng');
    print([figpath filesep char(strcat('fig_empirical_priors_', modSpace(m).name))], '-dsvg');
    close;

end

end

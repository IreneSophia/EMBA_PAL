function [] = plotEmpiricalPriors(modSpace, subDir, saveDir, midx)
% Plot empirical priors
% Based on hessetal_spirl_analysis/job_runner_paper_figs.m, lines
% 404 - 579 -> produces figure 5A in their paper
%
% !! CAUTION !!: This is HARDCODED for the models used in the EMBA project. 
%
% INPUT
%   modSpace            struct              model space
%
%   subDir              char array          subfolder containing estimates
%
%   saveDir             char array          base output directory
%
%   midx                double vector       model indices, [] for all
%-----------------------------------------------------------------------------
%
% Copyright (C) 2024 Anna Yurova, Irene Sophia Plank, LMU University Hospital;
%               based on the hessetal_spirl_analysis toolbox by Alex Hess (2024), TNU, ETHZ
%
%
% This file is part of the EMBA_HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. For further details, see the file LICENSE or <http://www.gnu.org/licenses/>.
%-----------------------------------------------------------------------------

empiricalPriors = load(fullfile(saveDir, subDir, 'params_priors'));
model = empiricalPriors.mod;

if isempty(midx)
    midx = 1:size(modSpace, 2);
end

%% FIG: empirical prior distributions

figpath = fullfile(saveDir, 'figures', subDir);
if ~exist(figpath, 'dir')
    mkdir(figpath);
end

% x grid
x_min = -50;
x_max = 50;
x = x_min:0.001:x_max;

for m = midx

    % reset the values
    prc = [];
    obs = [];

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
    count = 1;
    subplot(3,3,count)
    j = 1;
    plot(x, prc(j).y, 'k', 'LineWidth', 1)
    hold on
    plot(x, prc(j).y_prior, 'k--')
    % add participant data
    scatter(model(m).prc_est(:,j), zeros(length(model(m).prc_est(:,j)),1), ...
        20, 'o','MarkerEdgeColor','b','MarkerEdgeAlpha', 0.5, 'LineWidth', 1)
    hold off
    % figure out y limits
    upp = round(max([prc(j).y prc(j).y_prior])*1.05, 2);
    low = 0 - upp*0.05;
    ylim([low upp])
    % set x limits
    xlim([-15 5])
    if contains(modSpace(m).prc, "hgf")
        t = '$\omega_2$';
    else
        t = '$logit\alpha$';
    end
    title(t, 'interpreter','latex', 'FontSize', 20);

    % add a legend
    legend('empirical priors', 'default priors',...
        'Interpreter','latex', 'Position', [0.05 0.005 0.9 0.05], 'Orientation', 'horizontal', 'FontSize', 14)
    
    % plot om3
    if contains(modSpace(m).prc, "hgf")
        if modSpace(m).prc_config.n_levels > 2
            count = count + 1;
            subplot(3,3,count);
            j = 2;
            plot(x, prc(j).y, 'k', 'LineWidth', 1)
            hold on
            plot(x, prc(j).y_prior, 'k--')
            % add participant data
            scatter(model(m).prc_est(:,j), zeros(length(model(m).prc_est(:,j)),1), ...
                20, 'o','MarkerEdgeColor','b','MarkerEdgeAlpha', 0.5, 'LineWidth', 1)
            hold off
            % figure out y limits
            upp = round(max([prc(j).y prc(j).y_prior])*1.05, 2);
            low = 0 - upp*0.05;
            ylim([low upp])
            % set x limits
            xlim([-15 5])
            t = '$\omega_3$';
            title(t, 'interpreter','latex', 'FontSize', 20);
        end
    end

    %% RESPONSE MODEL

    for k = 1:(length(obs)-1)
        count = count + 1;
    
        % plot betas
        subplot(3,3,count)
        idx = find(round(obs(k).y, 3) > 0);
        yvalue = obs(k).y(idx);
        xvalue = x(idx);
        pvalue = obs(k).y_prior(idx);
        plot(xvalue, yvalue, 'k', 'LineWidth', 1)
        hold on
        plot(xvalue, pvalue, 'k--')
        % add participant data
        scatter(model(m).obs_est(:,k), zeros(length(model(m).obs_est(:,k)),1), ...
            20, 'o','MarkerEdgeColor','b','MarkerEdgeAlpha', 0.5, 'LineWidth', 1)
        hold off
        % figure out y limits
        upp = round(max([obs(k).y obs(k).y_prior])*1.05, 2);
        low = 0 - upp*0.05;
        ylim([low upp])
        t = sprintf('$\\beta_%d$', k-1);
        title(t, 'interpreter','latex', 'FontSize', 20);

    end
    
    % Plot Sigma (which is always the last parameter)
    count = count + 1;
    subplot(3,3,count)
    plot(x, obs(end).y, 'k', 'LineWidth', 1)
    hold on
    plot(x, obs(end).y_prior, 'k--')
    % add participant data
    scatter(model(m).obs_est(:,end), zeros(length(model(m).obs_est(:,end)),1), ...
        20, 'o','MarkerEdgeColor','b','MarkerEdgeAlpha', 0.5, 'LineWidth', 1)
    hold off
    % figure out y limits
    upp = round(max([obs(k).y obs(k).y_prior])*1.05, 2);
    low = 0 - upp*0.05;
    ylim([low upp])
    % set x limits
    xlim([-5 4])
    t = '$log(\Sigma)$';
    title(t, 'interpreter','latex', 'FontSize', 20);
    
    print([figpath filesep char(strcat('fig_empirical_priors_', modSpace(m).name))], '-dpng');
    print([figpath filesep char(strcat('fig_empirical_priors_', modSpace(m).name))], '-dsvg');
    close;

end

end

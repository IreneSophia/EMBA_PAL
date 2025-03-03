function [] = plotInitialPriorsPredDens(modSpace, u, saveDir)
% Plot initial prior predictive densities:
% 'Initial_priors_ehgf_sample_om2_om3',
% 'Initial_priors_ehgf_prior_pred_dens'
% Based on hessetal_spirl_analysis/utils/spirl_init_priors.m,lines 243-345
%
% !! CAUTION !!: This is HARDCODED for the models used in the EMBA project. 
%
% INPUT
%   modSpace            struct              model space
%                                           
%   u                   matrix              inputs at every trial
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

options = setupOpt_and_randSeed();

for m = 1:size(modSpace, 2)

    %% initial priors
    c.c_prc = modSpace(m).prc_config;
    prc_vect_nat = c.c_prc.transp_prc_fun(c, c.c_prc.priormus);
    prc_idx = modSpace(m).prc_idx;
    
    c.c_obs = eval(char(strcat(modSpace(m).obs, '_config')));
    obs_vect_nat = c.c_obs.transp_obs_fun(c, c.c_obs.priormus);
    
    %% sample omega priors and create simulated data

    % create empty omega matrix
    om = NaN(100,length(modSpace(m).prc_idx));
    
    % define 100 parameter values (1 priormu + 99 samples)
    n_samples = 100;
    
    % generate synthetic data based on prior mus and sas
    om(1,:) = c.c_prc.priormus(modSpace(m).prc_idx);
    omsas   = c.c_prc.priorsas(modSpace(m).prc_idx);
    prc_vect_nat(modSpace(m).prc_idx) = om(1,:);
    simData(1) = tapas_simModel(u,...
        char(modSpace(m).prc),...
        prc_vect_nat,...
        char(modSpace(m).obs),...
        obs_vect_nat,...
        options.rng.settings.State(1)); % seed

    % randomly sample the omegas based on their priors
    count   = 2;

    % do this until we have 100 synthetic data sets
    while count <= n_samples

        % sample new omegas
        for j = 1:length(om(1,:))
            om(count,j) = normrnd(om(1,j), sqrt(omsas(j)));
        end
        
        try

            % add to model
            prc_vect_nat(modSpace(m).prc_idx) = om(count,:);

            % generate synthetic data
            simData(count) = tapas_simModel(u,...
                char(modSpace(m).prc),...
                prc_vect_nat,...
                char(modSpace(m).obs),...
                obs_vect_nat); %,...
                %options.rng.settings.State()); %seed

            % if it worked increase the counter
            count = count + 1;
            
        catch

            warning('Unstable for omega: %.2f and %.2f', ...
                om(count,1), om(count,2))

        end
    
    end
    
    %% plot prior predictive densities of simulated data

    figure
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    
    % loop through the muhats
    for j = 1:modSpace(m).prc_config.n_levels
        x = [];
        subplot(3,2,7-2*j)
        % limits of plot
        if j == 1
            upper = 1;
            lower = 0;
            lim_upper = 1.1;
            lim_lower = -0.1;
        elseif j == 2
            upper = 10;
            lower = -10;
            lim_upper = 12;
            lim_lower = -12;
        elseif j == 3
            upper = 2;
            lower = -6;
            lim_upper = 2.8;
            lim_lower = -6.8;
        end
        % input
        input = u(:,1)*upper;
        input(input == 0) = lower;
        scatter(1:length(input), input, 5, 'o','MarkerEdgeColor','k','MarkerEdgeAlpha', 0.6, 'LineWidth', 0.1)
        hold on
        % plot the muhat traces
        for n = 1:n_samples
            x = [x simData(n).traj.muhat(:,j)];
            if n == 1
                plot(simData(n).traj.muhat(:,j), 'Color', [0 0 1 1], 'LineWidth', 2)
            else
                plot(simData(n).traj.muhat(:,j), 'Color', [0 0 1 0.2])
            end
        end
        % plot the median muhat
        plot(median(x, 2), 'Color', 'black', 'LineWidth', 3)
        xlabel('trials', 'FontSize', 20)
        ytxt = ['$\hat{\mu}_{', num2str(j), '}$'];
        ylabel(ytxt, 'Interpreter', 'Latex', 'FontSize', 20)
        txt = ['$\hat{\mu}_' num2str(j) '$'];
        title(txt, 'Interpreter', 'Latex', 'FontSize', 20)
        ylim([lim_lower lim_upper])
    end
    % loop through the sigmahats
    for k = 1:modSpace(m).prc_config.n_levels
        x = [];
        subplot(3,2,8-2*k)
        % limits of plot
        lower = 0;
        if k == 1
            upper = 0.3;
            lim_upper = 0.4;
            lim_lower = -0.1;
        elseif k == 2
            upper = 200;
            lim_upper = 220;
            lim_lower = -20;
        elseif k == 3
            upper = 200;
            lim_upper = 220;
            lim_lower = -20;
        end
        % input
        input = u(:,1)*upper;
        scatter(1:length(input), input, 5, 'o','MarkerEdgeColor','k','MarkerEdgeAlpha', 0.6, 'LineWidth', 0.1)
        hold on
        % plot the sigmahat trace
        for n = 1:n_samples
            x = [x simData(n).traj.sahat(:,k)];
            if n == 1
                plot(simData(n).traj.sahat(:,k), 'Color', [0 0 1 1], 'LineWidth', 2)
            else
                plot(simData(n).traj.sahat(:,k), 'Color', [0 0 1 0.2])
            end
        end
        % plot the median sigmas
        plot(median(x, 2), 'Color', 'black', 'LineWidth', 3)
        xlabel('trials')
        ytxt = ['$\hat{\sigma}_{', num2str(k), '}$'];
        ylabel(ytxt, 'Interpreter', 'Latex', 'FontSize', 20)
        txt = ['$\hat{\sigma}_' num2str(k) '$'];
        title(txt, 'Interpreter', 'Latex', 'FontSize', 20)
        ylim([lim_lower lim_upper])
    end

    figdir = fullfile(saveDir, 'figures', 'priors');
    if ~exist(figdir, 'dir')
        mkdir(figdir)
    end
    print([figdir filesep 'Initial_priors_ehgf_prior_pred_dens_' convertStringsToChars(modSpace(m).name)], '-dsvg');
    print([figdir filesep 'Initial_priors_ehgf_prior_pred_dens_' convertStringsToChars(modSpace(m).name)], '-dpng');
    close;

end

end

function [] = plotLogRTsModel(modSpace, subDir, saveDir, nSubs, midx)
% Plot logRT fits: 'logRT_model'
% Based on hessetal_spirl_analysis/job_runner_plot_results.m, lines
% 311-405
%
% !! CAUTION !!: This is HARDCODED for the models used in the EMBA project. 
%
%
% INPUT
%   modSpace            struct              model space
%
%   subDir              char array          subdirectory containing results
%
%   saveDir             char array          base output directory
%
%   nSubs               double              no of subs to plot, 0 for all
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

nModels = size(modSpace, 2);

if isempty(midx)
    midx = 1:nModels;
end

if strcmp(subDir, 'pilots') || strcmp(subDir, 'main')
    load(fullfile(saveDir, subDir, 'full_results'), 'res');
elseif strcmp(subDir, 'sim')
    load(fullfile(saveDir, subDir, 'recovery_analysis'), 'est');
    res = struct();
    for m = 1:size(est,1)
        for n = 1:size(est,2)
            res.est(m,n) = est(m,n,m).data;
        end
    end
    subDir = 'sim_rec';
end
nTrials = numel(res.est(1,1).y);
if nSubs == 0
    nSubs   = size(res.est, 2);
end
 

%% PILOT: plot rt trajectories & fits (single-sub)

nr_width  = 1;
nr_height = 4;
nr_total  = nr_width*nr_height;
nr_plots  = ceil(nSubs/nr_total);
for m = midx
    for p = 1:nr_plots
        figure
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        for i = 1:nr_total
            sub = i + (p-1)*nr_total;
            if sub > nSubs
                continue
            end
            subplot(nr_height,nr_width,i)
            % figure out subject-specific ylims
            lim_lower = floor((min(res.est(m,sub).y)-0.2)*2)/2;
            lim_upper = ceil((max(res.est(m,sub).y)+0.2)*2)/2;
            % adjust input to them
            input   = res.est(m,sub).u(:,1);
            input(input == 1) = lim_upper - 0.2;
            input(input == 0) = lim_lower + 0.2;
            % plot the volatile phase
            fill([73 264 264 73], [lim_lower lim_lower lim_upper lim_upper], ...
                     'b', 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
            hold on
            % plot the prevolatile phase
            fill([1 73 73 1], [lim_lower lim_lower lim_upper lim_upper], ...
                     'r', 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
            % plot the prevolatile phase
            fill([264 nTrials nTrials 264], [lim_lower lim_lower lim_upper lim_upper], ...
                     'r', 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
            % plot the fitted and real log rts
            plot(res.est(m,sub).y, 'Linewidth', 2); %, 'color', options.col.tnub)
            plot(res.est(m,sub).optim.yhat, 'r', 'LineWidth',2);
            scatter(1:length(input), input, 5, 'o','MarkerEdgeColor','k','MarkerEdgeAlpha', 0.6, 'LineWidth', 0.1)
            xlim([0 nTrials])
            ylim([lim_lower lim_upper])
            ylabel('logRT [ms]', 'FontSize', 14)
            xlabel('trials', 'FontSize', 14)
            % add the lme of the fit to the title
            txt = sprintf('%s sub %d (LME: %.03f)', subDir, sub, res.est(m,sub).optim.LME);
            title(txt, 'FontSize', 20)
        end
        legend('', '', '', '$log(y_{rt})$', '$log(\hat{y}_{rt})$', 'Interpreter','latex', 'Position', [0.94 0.48 0.03 0.07], 'FontSize', 20)
        figdir = fullfile(saveDir, 'figures', [subDir '_logRTs']); 
        if ~exist(figdir, 'dir')
            mkdir(figdir)
        end
        print(strcat(figdir, filesep, 'logRT_model', modSpace(m).name, '_', num2str(p)), '-dpng');
        print(strcat(figdir, filesep, 'logRT_model', modSpace(m).name, '_', num2str(p)), '-dsvg');
        close;
    end
end

end

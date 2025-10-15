function [] = plotSimData(modSpace, nSim, saveDir)
% Plot empirical priors: 'fig5C_empirical_prior_pred_dens_log_yhat_rt_M'
% Based on hessetal_spirl_analysis/job_runner_paper_figs.m, lines
% 621-653
% 
% !! CAUTION !!: This is HARDCODED for the models used in the EMBA project. 
%
% INPUT
%   nModels             integer             number of models
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

nMod = size(modSpace, 2);
sim = load(fullfile(saveDir, 'sim', 'simulated_data'));
options = setupOpt_and_randSeed();

%% get number of trials

nTrials = length(sim.sub(1).data.u(:,1));

%% specify colors
options.col.idx = ones(1,nTrials);

%% get screensize & hardcode figure sizes
scrsz = get(0,'screenSize');
% hard code figure sizes
outerpos1 = [scrsz(3),0.4*scrsz(4),700,250];

% x grid
x_min = -30;
x_max = 30;
x = x_min:0.1:x_max;
 
%% FIG: prior pred densities under empirical priors, log(yhat_rt) (M1)
% create figure
figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 1]);
% collect simulated logRTs
for m=1:nMod
    logRT_sim = NaN(nTrials,nSim);
    for n = 1:nSim
        logRT_sim(:,n) = sim.sub(n,m).data.y;
    end
    % ignore infinite values when computing the median
    logRT_sim(isinf(logRT_sim)) = nan;
    avg_logRT_sim = median(logRT_sim, 2, 'omitnan');
    subplot(4,2,m)
    plot(logRT_sim, 'color', [0 0 1 0.1])
    hold on
    % Time window is HARDCODED HERE!
    plot(1:nTrials, ones(nTrials,1)*log(1700), '--k')
    plot(1:nTrials, ones(nTrials,1)*log(100), '--k')
    plot(avg_logRT_sim, 'color', [0 0 1], 'LineWidth', 3)
    ylim([3 10])
    ax = gca; ax.FontSize = 14;
    xlim([1 nTrials]);
    ylabel('$log(y_{RT})$', 'Interpreter','Latex')
    xlabel('trials')
    title(strcat(modSpace(m).name, ' pilot simulation results'));
    hold off

end

figdir = fullfile(saveDir, 'figures', 'sim_data');
if ~exist(figdir, 'dir')
   mkdir(figdir)
end
print(fullfile(figdir, 'logRTs'), '-dpng');
print(fullfile(figdir, 'logRTs'), '-dsvg');
close;

end

function [] = plotTrajForMain(modSpace, m, parallel, local_cores, saveDir)
% 'binary_model' plots (or trajectory plots) for all simulated subjects
% Based on hessetal_spirl_analysis/job_runner_plot_results.m, lines
% 919-926
%
% INPUT
%   modSpace            struct              model space
%
%   parallel            0 or 1              parallelizaion switch:
%                                           1 - parallel
%                                           0 - serial
%
%   local_cores         double              how many cores to use for parallel
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

%% REC: plot continuous trajectories of the models per subject

load(fullfile(saveDir, 'main', 'full_results.mat'), 'res');

% Optional plotting of standard deviation for upper levels (true or false)
plotsd = true;

if (parallel)
    % request cores
    local = parcluster('local'); local.NumWorkers = local_cores;
    pool = parpool(local, local.NumWorkers);
    maxNumWorkers = local_cores;
else
    maxNumWorkers = 0;
end

parfor (sub=1:size(res.est,2), maxNumWorkers)
    emba_ehgf_binary_pu_tbt_plotTraj(res.est(m,sub), plotsd)
    hold on
    annotation('textbox', [0.48 0.88 0.8 0.1], ...
        'String', sprintf('LME: %.3f', res.est(m,sub).optim.LME), ...
        'Color', [1 0.5 0], ...
        'FontWeight', 'bold', ...
        'EdgeColor', 'none')
    hold off
    figdir = fullfile(saveDir, 'figures', 'main_traj');
    if ~exist(figdir, 'dir')
        mkdir(figdir)
    end
    print(strcat(figdir, filesep, 'traj', '_sub-', num2str(res.est(m,sub).subID), '_', modSpace(m).name), '-dpng');
    print(strcat(figdir, filesep, 'traj', '_sub-', num2str(res.est(m,sub).subID), '_', modSpace(m).name), '-dsvg');
    close;
end

% shuting down cores
if (parallel)
    pool.delete()
end


end

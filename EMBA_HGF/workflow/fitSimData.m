function [] = fitSimData(modSpaceRobMean, nSim, parallel, local_cores, saveDir)
% Fit all models for each simulated subject.
% Based on hessetal_spirl_analysis/job_runner_5_sim_data_modinv.m
%
% INPUT
%   modSpaceRobMean     struct              updated model space with
%                                           empirical priors
%
%   nSim                integer             number of simulated subjects
%
%   parallel            0 or 1              parallelizaion switch:
%                                           1 - parallel
%                                           0 - serial
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

% Load simulated data
sim = load(fullfile(saveDir, 'sim', 'simulated_data'));

% Number of models
nModels = size(modSpaceRobMean, 2);

if (parallel)
    % request cores
    local = parcluster('local'); local.NumWorkers = local_cores;
    pool = parpool(local, local.NumWorkers);
    maxNumWorkers = local_cores;
else
    maxNumWorkers = 0;
end

% Fit simulated data
parfor (n=1:nSim, maxNumWorkers)
    for m=1:nModels
        for i=1:nModels
            invertModelSim(n, m, i, modSpaceRobMean, sim, saveDir)
        end
    end
end

if (parallel)
    % shuting down cores
    pool.delete()
end

end

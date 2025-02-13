function [] = fitPilotData(modSpace, nPilots, pilot, parallel, local_cores, saveDir)
% Fit all models for each pilot subject.
% Based on hessetal_spirl_analysis/job_runner_2_pilots_modinv.m
%
% INPUT
%   modSpace            struct              model space
%
%   nPilots             integer             number of pilot participants
%
%   pilot               struct              pilot data
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

disp('Fitting pilot data')

% Number of models
nModels = size(modSpace, 2);

if (parallel)
    % request cores
    local = parcluster('local'); local.NumWorkers = local_cores;
    pool = parpool(local, local.NumWorkers);
    maxNumWorkers = local_cores;
else
    maxNumWorkers = 0;
end

parfor (n=1:nPilots, maxNumWorkers)
    for m = 1:nModels
        fprintf('------------------ %d ------------------', m)
        invertModelPilot(n, m, modSpace, pilot, saveDir);
    end
end

% Accumulate the results. Necessary for some of the plotting. Not the most
% efficient way to do it. Can be commented out if the corresponding plots
% are not used.
accumulateForPlotting('pilots', modSpace, saveDir)

% shuting down cores
if (parallel)
    pool.delete()
end

disp('Pilot data fitting done.')
end

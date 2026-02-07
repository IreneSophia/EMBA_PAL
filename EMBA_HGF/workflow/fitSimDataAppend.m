function [] = fitSimDataAppend(modSpaceAll, nSim, parallel, local_cores, saveDir, revDir)
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
%   local_cores         double              how many cores to use for parallel
%
%   saveDir             char array          base output directory
%
%   revDir              char array          new output directory
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
sim = load(fullfile(revDir, 'sim', 'simulated_data'));

% Number of models
nModels = size(modSpaceAll, 2);

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
            % check if the model already exists in saveDir
            save_path = fullfile(saveDir, 'sim', ['sim_sub', num2str(n)]);
            filename = [save_path filesep 'sim_mod_', convertStringsToChars(modSpaceAll(m).name),...
                    '_est_mod_', convertStringsToChars(modSpaceAll(i).name) '.mat'];
            if exist(filename, 'file')
                % copy it to the new directory
                dirnew  = fullfile(revDir, 'sim', ['sim_sub', num2str(n)]);
                if ~exist(dirnew, 'dir')
                   mkdir(dirnew)
                end
                filenew = [dirnew filesep 'sim_mod_', convertStringsToChars(modSpaceAll(m).name),...
                    '_est_mod_', convertStringsToChars(modSpaceAll(i).name) '.mat'];
                copyfile(filename, filenew);

            else 
                % run it and save it in the new directory
                invertModelSim(n, m, i, modSpaceAll, sim, revDir)
            end
        end
    end
end

if (parallel)
    % shuting down cores
    pool.delete()
end

end

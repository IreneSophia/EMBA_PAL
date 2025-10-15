function [model] = computeEmpiricalPriors(modSpace, nSubs, subDir, saveDir)

% Compute empirical priors from based on the pilot data fitting.
% Based on hessetal_spirl_analysis/job_runner_3_pilots_priors.m
%
% INPUT
%   modSpace            struct              model space
%
%   nSubs               integer             number of participants
%
%   subDir              char array          subfolder containing estimates
%
%   saveDir             char array          base output directory
%
%-----------------------------------------------------------------------------
%
% Copyright (C) 2024 Anna Yurova, Irene Sophia Plank, LMU University Hospital;
%               based on the hessetal_spirl_analysis toolbox by Alex Hess (2024), TNU, ETHZ
%
%
% This file is part of the EMBA_HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. For further details, see the file LICENSE or <http://www.gnu.org/licenses/>.
%-----------------------------------------------------------------------------

% Structure to store intermediate data
model = struct();

% Get a list of all the subjects in this folder
ls_subs = dir((fullfile(saveDir, subDir, 'sub*')));

% Aggregate estimated parameter values of all subjects for each model 
for n = 1:nSubs
    for m = 1:size(modSpace, 2)
        fprintf('current iteration: n=%1.0f, m=%1.0f \n', n,m);
        est = load(fullfile(ls_subs(n).folder, ls_subs(n).name, ...
            ['est_mod_', convertStringsToChars(modSpace(m).name)]));
        
        % prc model
        for j = 1:size(modSpace(m).prc_idx,2)

            % Save the estimated parameter value
            model(m).prc_est(n,j) = est.p_prc.ptrans(modSpace(m).prc_idx(j));

            % Save initial priors (they are the same for all subjects)
            if n == 1
                model(m).prc_priormus(j) = modSpace(m).prc_config.priormus(modSpace(m).prc_idx(j));
                model(m).prc_priorsas(j) = modSpace(m).prc_config.priorsas(modSpace(m).prc_idx(j));
            end
        end
        
        % obs model
        for k = 1:size(modSpace(m).obs_idx,2)

            % Save the estimated parameter value
            model(m).obs_est(n,k) = est.p_obs.ptrans(modSpace(m).obs_idx(k));

            % Save initial priors (they are the same for all subjects)
            if n == 1
                model(m).obs_priormus(k) = modSpace(m).obs_config.priormus(modSpace(m).obs_idx(k));
                model(m).obs_priorsas(k) = modSpace(m).obs_config.priorsas(modSpace(m).obs_idx(k));
            end
        end           
    end
end

%% estimate pilot priors

if strcmp(subDir, 'pilots')
    for m = 1:size(modSpace, 2)
        % prc
        model(m).prc_robmean = NaN(size(model(m).prc_priormus));
        model(m).prc_robvar  = NaN(size(model(m).prc_priorsas));
    
        % Compute robust variance and mean from for each free parameter based
        % on aggregated parameter values from all subjects
        for j = 1:size(modSpace(m).prc_idx,2)
            [var, model(m).prc_robmean(1,j)] = robustcov(model(m).prc_est(:,j));
            model(m).prc_robvar(1,j) = var;
        end
        
        % obs
        model(m).obs_robmean = NaN(size(model(m).obs_priormus));
        model(m).obs_robvar  = NaN(size(model(m).obs_priorsas));
    
        % Compute robust variance and mean from for each free parameter based
        % on aggregated parameter values from all subjects
        for k = 1:size(modSpace(m).obs_idx,2)
            [var, model(m).obs_robmean(1,k)] = robustcov(model(m).obs_est(:,k));
            model(m).obs_robvar(1,k) = var;
        end
    end
else
    empiricalPriors = load(fullfile(saveDir, 'pilots', 'params_priors'));
    for m = 1:size(modSpace,2)
        model(m).prc_robmean = empiricalPriors.mod(m).prc_robmean;
        model(m).prc_robvar  = empiricalPriors.mod(m).prc_robvar;
        model(m).obs_robmean = empiricalPriors.mod(m).obs_robmean;
        model(m).obs_robvar  = empiricalPriors.mod(m).obs_robvar;
    end
end

% save to struct
empiricalPriors.mod = model;
save_path = fullfile(saveDir, subDir, ['params_priors']);
if ~exist(fullfile(saveDir, subDir), 'dir')
   mkdir(fullfile(saveDir, subDir))
end
save(save_path, '-struct', 'empiricalPriors');

end

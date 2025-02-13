function [] = checkParameterRecovery(modSpaceRobMean, nSim, saveDir)
% Check if the fitting recovered parameters with which the simulated
% outputs were generated
% Based on hessetal_spirl_analysis/job_runner_6_rec_analysis.m
%
% INPUT
%   modSpaceRobMean     struct              updated model space with
%                                           empirical priors
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

% Load simulated data
sim = load(fullfile(saveDir, 'sim', 'simulated_data'));

% Number of models
nModels = size(modSpaceRobMean, 2);

% Create res struct
res = struct();
res.sim = sim;

%% Load and store results from model inversion
rec = struct();
for n = 1:nSim
    for m = 1:nModels
        fprintf('current iteration: n=%1.0f, m=%1.0f \n', n,m);

        for i = 1:nModels
            % Load results from model inversion
            rec.est(m,n,i).data = load(fullfile(saveDir, ...
                'sim', ['sim_sub', num2str(n)],...
                ['sim_mod_', convertStringsToChars(modSpaceRobMean(m).name),...
                '_est_mod_', convertStringsToChars(modSpaceRobMean(i).name)]));

            % Store LME in matrix
            rec.model(m).LME(n,i) = rec.est(m,n,i).data.optim.LME;
        end
        % Store parameter values (sim & est)
        rec.param.prc(m).sim(n,:) = res.sim.sub(n,m).input.prc.transInp(modSpaceRobMean(m).prc_idx);
        rec.param.obs(m).sim(n,:) = res.sim.sub(n,m).input.obs.transInp(modSpaceRobMean(m).obs_idx);
        rec.param.prc(m).est(n,:) = rec.est(m,n,m).data.p_prc.ptrans(modSpaceRobMean(m).prc_idx);
        rec.param.obs(m).est(n,:) = rec.est(m,n,m).data.p_obs.ptrans(modSpaceRobMean(m).obs_idx);
    end
end

%% Parameter recovery (Pearson's correlation coefficient)
% Compute correlation between the parameters used for the simulation
% and the ones recovered after the fitting
for m = 1:nModels

    % Perception model
    [prc_coef, prc_p] = corr(rec.param.prc(m).sim, rec.param.prc(m).est);
    rec.param.prc(m).pcc = diag(prc_coef);
    rec.param.prc(m).pval = diag(prc_p);

    % Observation model
    [obs_coef, obs_p] = corr(rec.param.obs(m).sim, rec.param.obs(m).est);
    rec.param.obs(m).pcc = diag(obs_coef);
    rec.param.obs(m).pval = diag(obs_p);
end

%% Save results as struct
save_path = fullfile(saveDir, 'sim');
if ~exist(save_path, 'dir')
   mkdir(save_path)
end
save([save_path filesep 'recovery_analysis'], '-struct', 'rec');

end

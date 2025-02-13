function [modSpaceRobMean] = updateModSpace(modSpace, saveDir)
% Update model space to include empirical priors.
% Based on hessetal_spirl_analysis/job_runner_3_pilots_priors.m
%
% INPUT
%   modSpace            struct              model space
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

% Load empirical priors
empiricalPriors = load(fullfile(saveDir, 'pilots', ['pilot_priors']));
model = empiricalPriors.mod;

% Initialize model space with the old one
modSpaceRobMean = modSpace;

% Set pilot priors to the empirical priors, rounded up to 4th digit
for m = 1:size(modSpace, 2)
    % prc
    for j = 1:size(modSpace(m).prc_idx,2)
        modSpaceRobMean(m).prc_config.priormus(modSpace(m).prc_idx(j)) = ...
            round(model(m).prc_robmean(j), 4);
        modSpaceRobMean(m).prc_config.priorsas(modSpace(m).prc_idx(j)) = ...
            round(model(m).prc_robvar(j), 4);
    end

    % Align parameter fields with the explicit prior definitions with the 
    % content of the vectors c.priormus and c.priorsas 
    modSpaceRobMean(m).prc_config = tapas_align_priors_fields(...
        modSpaceRobMean(m).prc_config);
    % obs
    for k = 1:size(modSpace(m).obs_idx,2)
        modSpaceRobMean(m).obs_config.priormus(modSpace(m).obs_idx(k)) = ...
            round(model(m).obs_robmean(k), 4);
        modSpaceRobMean(m).obs_config.priorsas(modSpace(m).obs_idx(k)) = ...
            round(model(m).obs_robvar(k), 4);
    end

    % Align parameter fields with the explicit prior definitions with the 
    % content of the vectors c.priormus and c.priorsas
    modSpaceRobMean(m).obs_config = tapas_align_priors_fields(...
        modSpaceRobMean(m).obs_config);
end
end

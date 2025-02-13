function [ModelSpace] = setupModelSpace(obsNamesDisp, prcNamesDisp, obsNamesList, prcNamesList)
% Sets up the initial model space. 
% Based on hessetal_spirl_analysis/utils/spirl_setup_model_space.m
%
% INPUT
%   obsNamesDisp        array of strings    array of display names for obs
% 
%   prcNamesDisp        array of strings    array of display names for prc
%
%   obsNamesList        array of strings    array of observation model
%                                           names
%
%   prcNamesList        array of strings    array of perception model
%                                           names
%-----------------------------------------------------------------------------
% 
% Copyright (C) 2024 Anna Yurova, Irene Sophia Plank, LMU University Hospital; 
%               based on the hessetal_spirl_analysis toolbox by Alex Hess (2024), TNU, ETHZ
%
%
% This file is part of the EMBA_HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% LICENSE or <http://www.gnu.org/licenses/>.


%% init struct
ModelSpace = struct();

% counter variable
x = 0;

% loop through the perception models
for i = 1:numel(prcNamesList)

    % loop through the observation models
    for j = 1:numel(obsNamesList)

        % increase counter 
        x = x + 1;

        % Combine model names
        ModelSpace(x).name = prcNamesDisp(i) + "-" + obsNamesDisp(j);

        % Setup the perception model
        ModelSpace(x).prc = prcNamesList(i);
        ModelSpace(x).prc_config = eval(strcat(prcNamesList(i), '_config'));
    
        % Setup the observation model
        ModelSpace(x).obs = obsNamesList(j);
        ModelSpace(x).obs_config = eval(strcat(obsNamesList(j), '_config'));

    end

end

%% Find free parameters
% Note that it relies on the fact, that in the perception model, the
% variance for all parameters, but omegas, are set to zero!

% Loop through the models
for i = 1:size(ModelSpace,2)

    % Perception model

    % Get the vector prior variances of the perception model
    prc_idx = ModelSpace(i).prc_config.priorsas;

    % Set the values that are NaN to zero
    prc_idx(isnan(prc_idx)) = 0;

    % Find the non-zero indexes, which correspond to the free parameters
    ModelSpace(i).prc_idx = find(prc_idx);


    % Observation model

    % Get the vector prior variances of the observation model
    obs_idx = ModelSpace(i).obs_config.priorsas;

    % Set the values that are NaN to zero
    obs_idx(isnan(obs_idx)) = 0;

    % Find the non-zero indexes, which correspond to the free parameters
    ModelSpace(i).obs_idx = find(obs_idx);

end

end

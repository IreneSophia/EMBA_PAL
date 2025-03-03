function [lawModel, modSpace] = addLawModelSpace(modSpace)
% Add the model from Lawson et al. (2017) to the model space
%
% INPUT
% 
%   modSpace            struct              model space
% 
%-----------------------------------------------------------------------------
% 
% Copyright (C) 2024 Anna Yurova, Irene Sophia Plank, LMU University Hospital
%
%
% This file is part of the EMBA_HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% LICENSE or <http://www.gnu.org/licenses/>.

lawModel = size(modSpace, 2) + 1;
modSpace(lawModel).name = "Lawson2017";
modSpace(lawModel).obs = "tapas_logrt_linear_binary";
modSpace(lawModel).obs_config = eval(strcat(modSpace(end).obs, '_config'));
modSpace(lawModel).obs_idx = find(~isnan(modSpace(end).obs_config.priorsas) & modSpace(end).obs_config.priorsas ~= 0);
modSpace(lawModel).prc = "tapas_hgf_binary_pu_tbt";
modSpace(lawModel).prc_config = eval(strcat(modSpace(end).prc, '_config'));
modSpace(lawModel).prc_idx = find(~isnan(modSpace(end).prc_config.priorsas) & modSpace(end).prc_config.priorsas ~= 0);
% adjust priors of zeta:
modSpace(lawModel).obs_config.logzemu = 1;
modSpace(lawModel).obs_config.logzesa = 0.6;
modSpace(lawModel).obs_config.priormus(end) = modSpace(lawModel).obs_config.logzemu;
modSpace(lawModel).obs_config.priorsas(end) = modSpace(lawModel).obs_config.logzesa;

end

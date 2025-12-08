function [obsNamesDisp, prcNamesDisp, obsNamesList, prcNamesList] = modelNames()
% Specify model types in this file. 
% The workflow assumes that the perception model is the same for all
% observation models.
%
% OUTPUT
%   modelDisplayList    array of strings    array of model display names, which
%                                           will be appended to the output file names.
%
%   obsNamesList        array of strings    array of observation model
%                                           names
%
%   prcName             char array          perception model name

%-----------------------------------------------------------------------------
%
% Copyright (C) 2024 Anna Yurova, Irene Sophia Plank, LMU University Hospital;
%               based on the hessetal_spirl_analysis toolbox by Alex Hess (2024), TNU, ETHZ
%
%
% This file is part of the EMBA_HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. For further details, see the file LICENSE or <http://www.gnu.org/licenses/>.
%-----------------------------------------------------------------------------

% List the display names of the observation models.
% It is best to keep them short, as these strings will be appended to
% output file names, also better not to use spaces or special characters
obsNamesDisp = ["L17", "L21"];

% List the names of observation models
% It is necessary to use DOUBLE QUOTATION MARKS "..."!
% Make sure it is the same number as the number of the display names!
 obsNamesList = [...
     "tapas_logrt_linear_binary", "emba_logrt_linear_binary_C"];
if (numel(obsNamesDisp) ~= numel(obsNamesList))
    error('The number of model names does not match the number of config names. Check your input!')
end

% Perception model name
% The workflow assumes that the perception model is the same for all
% observation models > this model needs to be in the tapas HGF folder
prcNamesList = ["emba_ehgf_binary_pu_tbt", "emba_hgf_binary_pu_tbt"]; %;

% List the display names of the perception models.
% It is best to keep them short, as these strings will be appended to
% output file names, also better not to use spaces or special characters
prcNamesDisp = ["eHGF", "HGF"];
if (numel(prcNamesDisp) ~= numel(prcNamesList))
    error('The number of model names does not match the number of config names. Check your input!')
end

end

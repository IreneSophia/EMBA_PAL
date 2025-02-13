function [] = invertModelPilot(n, m, modSpace, pilot, saveDir)
% Fit model with the index m to the data of the subject n
% Based on hessetal_spirl_analysis/job_runner_5_sim_data_modinv.m
%
% INPUT
%   n                   integer             subject id
%
%   m                   integer             id of the model to be fitted
%
%   modSpace            struct              model space 
%
%   pilot               struct              pilot data
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

% Set the random number generator options and the optimization algorithm
% This function is called many times through out the code and resets the
% seed for the random number generator every time
options = setupOpt_and_randSeed();

% Invert the model
est = tapas_fitModel(pilot.data(n).y,... % responses
            pilot.data(n).u,... % input sequence
            modSpace(m).prc_config,... %Prc fitting model
            modSpace(m).obs_config,... %Obs fitting model
            options.opt_config); %opt algo

% Save model fit as struct
save_path = fullfile(saveDir, 'pilots', ['sub', num2str(n)],...
    ['est_mod_', convertStringsToChars(modSpace(m).name)]);
if ~exist(fullfile(saveDir, 'pilots', ['sub', num2str(n)]), 'dir')
   mkdir(fullfile(saveDir, 'pilots', ['sub', num2str(n)]))
end
save(save_path, '-struct', 'est');

end 

function [] = invertModel(sub, modSpace, saveDir)
% Fit model with the index m to the data of the subject n
% Based on hessetal_spirl_analysis/job_runner_5_sim_data_modinv.m
%
% INPUT
%
%   sub                 struct              participant data
%
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

% Set the random number generator options and the optimization algorithm
% This function is called many times through out the code and resets the
% seed for the random number generator every time
options = setupOpt_and_randSeed();

% Loop through the models

for m = 1:length(modSpace)

    % check whether tbt
    if contains(modSpace(m).prc, "tbt")
        u = sub.u;
    else
        u = sub.u(:,1);
    end

    % Invert the model
    est = tapas_fitModel(sub.y,... % responses
                u,... % input sequence
                modSpace(m).prc_config,... %Prc fitting model
                modSpace(m).obs_config,... %Obs fitting model
                options.opt_config); %opt algo

    % add group information to the structure
    est.group = sub.group;
    est.subID = sub.subID;
    
    % Save model fit as struct
    save_path = fullfile(saveDir, 'main', ['sub' num2str(sub.subID)], "est_mod_" + modSpace(m).name);
    if ~exist(fullfile(saveDir, 'main', ['sub' num2str(sub.subID)]), 'dir')
       mkdir(fullfile(saveDir, 'main', ['sub' num2str(sub.subID)]))
    end
    save(save_path, '-struct', 'est');

end

end 

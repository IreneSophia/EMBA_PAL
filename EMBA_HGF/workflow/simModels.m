function [] = simModels(modSpace, nSim, u, saveDir, verbose, maxRep)
% Generate simulated synthetic data.
% Based on hessetal_spirl_analysis/job_runner_4_sim_setup.m
%
% INPUT
%   modSpace            struct              model space
%
%   nSim                integer             number of simulated subjects
%
%   u                   matrix              inputs at every trial
%
%   saveDir             char array          base output directory
%
%   verbose             boolean             whether to print all values
%
%   maxRep              integer             how often to retry if instable
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


% Set the random number generator options and the optimization algorithm
% This function is called many times through out the code and resets the
% seed for the random number generator every time
options = setupOpt_and_randSeed();

%% create synthetic data

% pre-allocate variables
sim   = struct();
input = struct();

% loop over simulated subjects & model space

for m = 1:size(modSpace,2)

    % For testing: reset the seed and the random number generator
    % every time for more reproducible results for each
    % model.
    %
    % For the comparison of the two identical models, this is
    % necessary in order to make sure that both times the random numbers
    % initialize from the same seed.
    %
    % In the original code, the order of the loops was reversed (first over
    % n and then over m) and there was no re-initialization. This has to be
    % restored for the direct comparisson with Hess et al. results.
    % Otherwise, the results will not match.
    options = setupOpt_and_randSeed();

    for n = 1:nSim

        % Set the parameter values to the prior values
        input.prc.transInp = modSpace(m).prc_config.priormus;
        input.obs.transInp = modSpace(m).obs_config.priormus;

        % Sample the values of the free parameters of the perception model
        for j = 1:size(modSpace(m).prc_idx,2)
            input.prc.transInp(modSpace(m).prc_idx(j)) = ...
                normrnd(modSpace(m).prc_config.priormus(modSpace(m).prc_idx(j)),...
                abs(sqrt(modSpace(m).prc_config.priorsas(modSpace(m).prc_idx(j)))));
        end

        % Sample the values of the free parameters of the observation model
        for k = 1:size(modSpace(m).obs_idx,2)
            input.obs.transInp(modSpace(m).obs_idx(k)) = ...
                normrnd(modSpace(m).obs_config.priormus(modSpace(m).obs_idx(k)),...
                abs(sqrt(modSpace(m).obs_config.priorsas(modSpace(m).obs_idx(k)))));
        end

        % create simulation input vectors (native space)
        c.c_prc = modSpace(m).prc_config;
        input.prc.nativeInp = modSpace(m).prc_config.transp_prc_fun(c, input.prc.transInp);
        c.c_obs = modSpace(m).obs_config;
        input.obs.nativeInp = modSpace(m).obs_config.transp_obs_fun(c, input.obs.transInp);

        % simulate responses
        stable = 0;
        count  = 0;
        
        while stable == 0 && count < maxRep
            
            % error message
            me = [];

            % increase the counter
            count = count + 1;

            % catches all errors and does not show the appropriate message. 
            % instead, it tries again until maxRep is reached
            try
                sim.sub(n,m).data = tapas_simModel(u,...
                    char(modSpace(m).prc),...
                    input.prc.nativeInp,...
                    char(modSpace(m).obs),...
                    input.obs.nativeInp,...
                    options.rng.settings.State(options.rng.idx, 1));
                stable = 1;
            catch me
                if count < maxRep
                    fprintf('simulation failed for Model %s, number of sim %1.0f \n', convertStringsToChars(modSpace(m).name), n);
                else
                    error('CAUTION: final simulation failed for Model %s, number of sim %1.0f \n', convertStringsToChars(modSpace(m).name), n);
                end
                if verbose 
                    fprintf('Prc Param Values: \n');
                    disp(input.prc.nativeInp);
                    fprintf('Obs Param Values: \n');
                    disp(input.obs.nativeInp);
                end
                % re-sample ONLY prc param values
                for j = 1:size(modSpace(m).prc_idx,2)
                    input.prc.transInp(modSpace(m).prc_idx(j)) = ...
                        normrnd(modSpace(m).prc_config.priormus(modSpace(m).prc_idx(j)),...
                        abs(sqrt(modSpace(m).prc_config.priorsas(modSpace(m).prc_idx(j)))));
                end
                input.prc.nativeInp = modSpace(m).prc_config.transp_prc_fun(c, input.prc.transInp);
            end
        end

        % save sim input
        sim.sub(n,m).input = input;

        % Update the rng state idx
        options.rng.idx = options.rng.idx+1;
        if options.rng.idx == (length(options.rng.settings.State)+1)
            options.rng.idx = 1;
        end
    end
end

options.rng.idx = 1;

%% save simulated data
save_path = fullfile(saveDir, 'sim', ['simulated_data']);
if ~exist(fullfile(saveDir, 'sim'), 'dir')
   mkdir(fullfile(saveDir, 'sim'))
end
save(save_path, '-struct', 'sim');

disp('simulated data successfully created.')

end

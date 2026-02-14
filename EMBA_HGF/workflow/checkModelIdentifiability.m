function [] = checkModelIdentifiability(nSim, modSpace, saveDir)
% For each model check whether it is the best LME fit for the 
% data, simulated by it 
% Based on hessetal_spirl_analysis/job_runner_6_rec_analysis.m
%
% INPUT
%   nSim                integer             number of simulated subjects
%
%   nModels             integer             number of models
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

nModels = size(modSpace,2);

% Load simulated data
sim = load(fullfile(saveDir, 'sim', 'simulated_data'));

% Create res struct
res = struct();
res.sim = sim;

% Create directory for figures
figdir = fullfile(saveDir, 'figures', 'model_identifiability');
if ~exist(figdir, 'dir')
   mkdir(figdir)
end

%% Load and store results from model inversion

for n = 1:nSim
    for m = 1:nModels
        fprintf('current iteration: n=%1.0f, m=%1.0f \n', n,m);
        for i = 1:nModels
            % load results from model inversion
            rec.est(m,n,i).data = load(fullfile(saveDir, ...
                'sim', ['sim_sub', num2str(n)],...
                ['sim_mod_', convertStringsToChars(modSpace(m).name), ...
                '_est_mod_', convertStringsToChars(modSpace(i).name)]));
            % store LME in matrix
            rec.model(m).LME(n,i) = rec.est(m,n,i).data.optim.LME;
        end     
     end
end

%% Model identifiability (LME Winner classification)

% pre-allocate
class.LMEwinner = NaN(nModels, nModels);
class.percLMEwinner = NaN(size(class.LMEwinner));

% Calculate the winner frequency for each data generating model
for m = 1:nModels

    % Find the maximum LME value among all simulated subjects
    [class.max(m).val, class.max(m).idx] = max(rec.model(m).LME, [], 2);

    % Compute the amound of the LME winners for each model
    for i = 1:nModels
        class.LMEwinner(m,i) = sum(class.max(m).idx==i);
    end

    % Percentage of LME winners with respect to the overall number of
    % simualated subjects
    class.percLMEwinner(m,:) = class.LMEwinner(m,:)./nSim;

    % Accuracy
    class.acc(m) = class.percLMEwinner(m,m);
end

% Balanced accuraccy
class.balacc = mean(class.acc);

% Chance threshold (inv binomial distr) based on Hess et al.
class.chancethr = binoinv(0.9, nSim, 1/nModels) / nSim;

% Save to struct
rec.class = class;

%% model identifiability (RFX BMS)

% pre-allocate
bmc.rfx.Ef = NaN(nModels, nModels);
bmc.rfx.ep = NaN(size(bmc.rfx.Ef));
bmc.rfx.pxp = NaN(size(bmc.rfx.Ef));

% toolbox settings for BMS results display
bmc.opt.verbose = false;
bmc.opt.DisplayWin = true;

txt = {};
for m = 1:size(modSpace, 2)
    % add to the title
    txt = [txt, char(modSpace(m).name)];
end
bmc.opt.modelNames = txt;

% run BMS for the data generating models
for m = 1:nModels
    L = rec.model(m).LME';
    % throw out participants for which L is -Inf or Inf
    L = L(:,abs(sum(L)) ~= Inf);
    [post, out] = VBA_groupBMC(L, bmc.opt); 
    figpath = fullfile(figdir,...
            ['model_ident_bmc_' convertStringsToChars(modSpace(m).name)]);
    print(figpath, '-dpng');
    print(figpath, '-dsvg');
    close;
    bmc.out{m} = out;
    bmc.post{m} = post;
    bmc.rfx.Ef(m,:)  = out.Ef';
    bmc.rfx.ep(m,:)  = out.ep;
    bmc.rfx.pxp(m,:) = out.pxp;
end

% save to struct
rec.bmc = bmc;

%% Save results as struct
save_path = fullfile(saveDir, 'sim');
if ~exist(save_path, 'dir')
   mkdir(save_path)
end
save([save_path filesep 'model_identifiability'], '-struct', 'rec');

end

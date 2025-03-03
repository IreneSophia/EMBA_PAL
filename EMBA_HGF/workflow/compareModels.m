function [out] = compareModels(modSpaceRobMean, saveDir)
% Model comparison based on the models run on the real data
% Based on hessetal_spirl_analysis/job_runner_12_main_results.m
%
% INPUT
%
%   modSpace            structure           model space information
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

% initialise some stuff
nModels = size(modSpaceRobMean,2);
res_dir = dir([saveDir filesep 'main' filesep 'sub*']);
nSubs   = length(res_dir);
res     = struct();

% set the filepaths and create folders is necessary
save_path = fullfile(saveDir, 'main');
if ~exist(save_path, 'dir')
   mkdir(save_path)
end
if ~exist([saveDir filesep 'figures' filesep 'main'], 'dir')
   mkdir([saveDir filesep 'figures' filesep 'main'])
end

%% Load and store results from model inversion

txt = {};
for m = 1:nModels
    % add to the title
    txt = [txt, char(modSpaceRobMean(m).name)];
    for n = 1:nSubs
        fprintf('current iteration: n=%1.0f, m=%1.0f \n', n, m);
        % load results from model inversion
        res.est(m,n).data = load(fullfile(...
            res_dir(n).folder,...
            res_dir(n).name,...
            ['est_mod_', convertStringsToChars(modSpaceRobMean(m).name), '.mat']));
        % store LME 
        res.LME(n,m) = res.est(m,n).data.optim.LME;
        % store the group to be able to the information
        group{n} = res.est(m,n).data.group;
     end
end

% compare all four
options.modelNames = txt;
options.verbose    = false;
[res.indiv.posterior, res.indiv.out] = VBA_groupBMC(res.LME', options);
out = res.indiv.out;

% save the plot
print([saveDir filesep 'figures' filesep 'main' filesep 'model_comparison'], '-dpng');
print([saveDir filesep 'figures' filesep 'main' filesep 'model_comparison'], '-dsvg');
close;

% save results as struct
save([save_path filesep 'model_comparison'], '-struct', 'res');

end

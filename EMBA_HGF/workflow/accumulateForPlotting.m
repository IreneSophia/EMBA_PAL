function [] = accumulateForPlotting(subDir, modSpace, saveDir)
% Accumulate the results for plotting. This is necessary to match the
% format of the hessetal_spirl_analysis toolbox.
%
% INPUT
%   subDir              char array          subdirectory
%
%   nModels             integer             number of models
%
%   saveDir             char array          base output directory
%-----------------------------------------------------------------------------
%
% Copyright (C) 2024 Anna Yurova, Irene Sophia Plank, LMU University Hospital;
%               based on the hessetal_spirl_analysis toolbox by Alex Hess (2024), TNU, ETHZ
%
% This file is part of the EMBA_HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. For further details, see the file LICENSE or <http://www.gnu.org/licenses/>.
%-----------------------------------------------------------------------------

% Structure for saving accumulated fitting results. Needed for plotting
res = struct();
nModels = size(modSpace, 2);

% specify subfolder
saveDir = fullfile(saveDir, subDir);

% Accumulate the results. Not the most efficient way to do it, but was
% necessary for plotting.
ls_files = dir([saveDir filesep 'sub*']);
for n=1:length(ls_files)
    for m = 1:nModels
         res.res.est(m,n) = load(fullfile(saveDir, ls_files(n).name,...
            ['est_mod_', convertStringsToChars(modSpace(m).name)]));
    end
end

% Save full results
save_path = fullfile(saveDir, 'full_results');
if ~exist(fullfile(saveDir), 'dir')
   mkdir(fullfile(saveDir))
end
save(save_path, '-struct', 'res');
end 

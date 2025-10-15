function [] = simPostRT(modSpace, subDir, saveDir, midx, nSim)
% Simulate nSim data sets per subject and save as a CSV for further
% plotting
%
%
% INPUT
%   modSpace            struct              model space
%
%   subDir              char array          subdirectory containing results
%
%   saveDir             char array          base output directory
%
%   midx                double vector       model indices, [] for all
%
%   nSim                double              number of sims per subject
%-----------------------------------------------------------------------------
%
% Copyright (C) 2024 Irene Sophia Plank, LMU University Hospital;
%
%
% This file is part of the EMBA_HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. For further details, see the file LICENSE or <http://www.gnu.org/licenses/>.
%-----------------------------------------------------------------------------

nModels = size(modSpace, 2);

if isempty(midx)
    midx = 1:nModels;
end

load(fullfile(saveDir, subDir, 'full_results'), 'res');
nTrials = numel(res.est(1,1).y);

nSubs = length(res.est);

%% Simulate data

for m = midx

    % initialise columns
    subID = nan(nSubs*nTrials*nSim,1);
    trl   = repmat(1:nTrials, 1, nSubs*nSim).';
    yhat  = nan(nSubs*nTrials*nSim,1);
    sim   = nan(nSubs*nTrials*nSim,1);

    for i = 1:nSubs

        % add the subject info
        first = (nTrials*nSim*(i-1)+1);
        last  = (nTrials*nSim*i);
        subID(first:last) = res.est(m,i).subID;

        % simulate some data
        for j = 1:nSim
            test = tapas_simModel(res.est(m, i).u,...
                char(res.est(m, i).c_prc.prc_fun),...
                res.est(m, i).p_prc.p,...
                char(res.est(m, i).c_obs.obs_fun),...
                res.est(m, i).p_obs.p);
            yhat((first+(nTrials*(j-1))):(first+(nTrials*j)-1)) = test.y;
            sim((first+(nTrials*(j-1))):(first+(nTrials*j)-1)) = j;
        end

    end
    tbl = table(subID, sim, trl, yhat);
    writetable(tbl, fullfile(saveDir, subDir, sprintf("ppc_data_%s.csv", modSpace(m).name)));
end

end

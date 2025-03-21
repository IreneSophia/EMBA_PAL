function [] = extractHGFcsv(modSpaceRobMean, saveDir, midx)
% Extract subject-specific parameters from winning model.
%
% INPUT
%
%   saveDir             char array          base output directory
%
%   midx                double vector       model indices, [] for all
%
%-----------------------------------------------------------------------------
%
% Copyright (C) 2024 Anna Yurova, Irene Sophia Plank, LMU University Hospital
%
% Please cite the papers from CITATION.txt when using the EMBA_HGF Toolbox.
%
% This file is part of the EMBA_HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. For further details, see the file LICENSE or <http://www.gnu.org/licenses/>.
%-----------------------------------------------------------------------------

% load the results and the model space including the empirical priors
res_dir = [saveDir filesep 'main' filesep];
load([res_dir 'full_results.mat'], 'res')

nModels = size(modSpaceRobMean, 2);
if isempty(midx)
    midx = 1:nModels;
end

for m = midx

    % only keep the current model
    est = res.est(m,:);
    
    % create the tables
    tbl  = table();  % wide with subjects
    ltbl = table();  % long with trial by trial prediction errors
    
    % get the trial number
    n_trials  = size(est(end).traj.muhat, 1);
    
    % loop trough the subjects
    for i = 1:length(est)

        %% Parameter table

        % extract subject info
        tbl.subID(i) = est(i).subID;
        tbl.diagnosis{i} = est(i).group;
        % extract model assessment
        tbl.LME(i)   = est(i).optim.LME; % more is better
        tbl.BIC(i)   = est(i).optim.BIC; % less is better
        tbl.AIC(i)   = est(i).optim.AIC; % less is better
        % extract betas from obs model
        for field = fieldnames(est(i).p_obs)'
            if length(est(i).p_obs.(field{1})) == 1
                tbl.(field{1})(i) = est(i).p_obs.(field{1});
            end
        end
        % extract omegas from perceptual model
        tbl.om2(i)   = est(i).p_prc.om(2);
        tbl.om3(i)   = est(i).p_prc.om(3);
        % add the info on the priors that were used for the obs model
        for field = fieldnames(est(i).c_obs)'
            if length(est(i).c_obs.(field{1})) == 1 && ~contains(field{1}, 'fun')
                tbl.(field{1})(i) = est(i).c_obs.(field{1});
            end
        end
        % add the info on the priors that were used for the prc model
        tbl.om2mu(i) = est(i).c_prc.ommu(2);
        tbl.om2sa(i) = est(i).c_prc.omsa(2);
        tbl.om3mu(i) = est(i).c_prc.ommu(3);
        tbl.om3sa(i) = est(i).c_prc.omsa(3);

        %% long table with trajectories
        
        idx = (height(ltbl)+1):(height(ltbl)+n_trials);
        ltbl.trl(idx)       = [1:n_trials].';
        ltbl.subID(idx)     = repmat(est(i).subID, n_trials, 1);
        ltbl.diagnosis(idx) = repmat(convertCharsToStrings(est(i).group), n_trials, 1);
        ltbl.eps2(idx)      = est(i).traj.epsi(:,2);
        ltbl.eps3(idx)      = est(i).traj.epsi(:,3);
        ltbl.mu1(idx)       = est(i).traj.mu(:,1);
        ltbl.mu2(idx)       = est(i).traj.mu(:,2);
        ltbl.mu3(idx)       = est(i).traj.mu(:,3);
        ltbl.sa2(idx)       = est(i).traj.sa(:,2);
        ltbl.sa3(idx)       = est(i).traj.sa(:,3);
        ltbl.muhat1(idx)    = est(i).traj.muhat(:,1);
        ltbl.muhat2(idx)    = est(i).traj.muhat(:,2);
        ltbl.muhat3(idx)    = est(i).traj.muhat(:,3);
        ltbl.sahat1(idx)    = est(i).traj.sahat(:,1);
        ltbl.sahat2(idx)    = est(i).traj.sahat(:,2);
        ltbl.sahat3(idx)    = est(i).traj.sahat(:,3);
        ltbl.da1(idx)       = est(i).traj.da(:,1);
        ltbl.da2(idx)       = est(i).traj.da(:,2);
        ltbl.da3(idx)       = est(i).traj.da(:,3);
        ltbl.psi2(idx)      = est(i).traj.psi(:,2);
        ltbl.psi3(idx)      = est(i).traj.psi(:,3);
        ltbl.alpha1(idx)    = est(i).traj.wt(:,1);
        ltbl.alpha2(idx)    = est(i).traj.wt(:,2);
        ltbl.alpha3(idx)    = est(i).traj.wt(:,3);
    
    end
    
    writetable(tbl,  sprintf('%s%s_results.csv', res_dir, modSpaceRobMean(m).name ))
    writetable(ltbl, sprintf('%s%s_traj.csv',    res_dir, modSpaceRobMean(m).name ))

end

end

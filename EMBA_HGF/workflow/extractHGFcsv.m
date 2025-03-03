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
    
    % create a table
    tbl = table();
    
    % get the trial number
    n_trials  = size(est(end).traj.muhat, 1);

    % initialise alpha matrices
    alpha2    = nan(n_trials,length(est));
    alpha3    = nan(n_trials,length(est));
    
    % loop trough the subjects
    for i = 1:length(est)
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
        % compute the learning rate following Lawson 2017, first alpha2
        for t = 2:n_trials
            alpha2(t,i)  = ...
                (est(i).traj.muhat(t, 1) - est(i).traj.muhat(t-1, 1))/...
                 est(i).traj.da(t, 1);
        end
        % then, alpha3
        for t = 2:n_trials
            alpha3(t,i)  = ...
                (est(i).traj.mu(t, 3) - est(i).traj.mu(t-1, 3))/...
                 est(i).traj.da(t, 2);
        end
        % add the learning rate updates to the table
        tbl.alpha2_pre2vol(i)  = median(alpha2(1:72,i),    "omitnan") - ...
            median(alpha2(73:144,i),  "omitnan");
        tbl.alpha2_vol2post(i) = median(alpha2(193:264,i), "omitnan") - ...
            median(alpha2(265:336,i), "omitnan");
        tbl.alpha3_pre2vol(i)  = median(alpha3(1:72,i),    "omitnan") - ...
            median(alpha3(73:144,i),  "omitnan");
        tbl.alpha3_vol2post(i) = median(alpha3(193:264,i), "omitnan") - ...
            median(alpha3(265:336,i), "omitnan");
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
    
    end

    % create a long table with alpha values
    atbl = table();
    for i = 1:height(tbl)
        idx = (height(atbl)+1):(height(atbl)+n_trials);
        atbl.trl(idx)       = [1:n_trials].';
        atbl.subID(idx)     = repmat(tbl.subID(i), n_trials, 1);
        atbl.diagnosis(idx) = repmat(tbl.diagnosis(i), n_trials, 1);
        atbl.alpha2(idx)    = alpha2(:,i);
        atbl.alpha3(idx)    = alpha3(:,i);
    end
    
    writetable(tbl,  sprintf('%s%s_results.csv', res_dir, modSpaceRobMean(m).name ))
    writetable(atbl, sprintf('%s%s_alphas.csv',  res_dir, modSpaceRobMean(m).name ))

end

end

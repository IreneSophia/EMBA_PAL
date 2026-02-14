function [] = plotAggLogRTsModelPhases(modSpace, subDir, saveDir, input)
% Plot comparison of real data and predicted data for conditions
%
% !! CAUTION !!: This is HARDCODED for the models used in the EMBA project. 
%
%
% INPUT
% 
%   modSpace            struct              model space
%
%   subDir              char array          subdirectory containing results
%
%   saveDir             char array          base output directory
%
%   input               structure           subject-specific input for model 
%-----------------------------------------------------------------------------
%
% Copyright (C) 2024 Anna Yurova, Irene Sophia Plank, LMU University Hospital;
%               based on the hessetal_spirl_analysis toolbox by Alex Hess (2024), TNU, ETHZ
%
%
% This file is part of the EMBA_HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. For further details, see the file LICENSE or <http://www.gnu.org/licenses/>.
%-----------------------------------------------------------------------------

nModels = size(modSpace, 2);

load(fullfile(saveDir, subDir, 'full_results'), 'res');

nTrials = numel(res.est(1,1).y);
nSubs   = size(res.est, 2);

figdir = fullfile(saveDir, 'figures', subDir); 
if ~exist(figdir, 'dir')
    mkdir(figdir)
end

%% Load condition indices

% load stimulus order > expectedness same in pilot and our data
tbl = readtable("PAL_scheme.csv");
coi = {'difficulty', 'expected'}; % conditions of interest

%% aggregate real and predicted data per partipant
phases = {'prevolatile', 'volatile', 'postvolatile'};

for m = 1:nModels
    % aggregation
    fout = cell(3,1);
    for p = 1:3
        phase = phases{p};
        out = struct();
        for field = coi
            opts      = unique(tbl.(field{1}));
            temp      = table();
            nrow      = nSubs*length(opts)*2;
            subs      = nan(nrow,1);
            cond      = cell(nrow,1);
            what      = subs;
            data      = subs;
            count     = 1;
            for i = 1:nSubs
                for j = 1:length(opts)
                    % find the corresponding indices with difficulty either in
                    % the u itself > different in Lawson and our study < or in
                    % the planned scheme for the models without difficulty
                    if strcmp(field{1}, 'difficulty')
                        if strcmp(opts{j}, 'easy')
                            idx = intersect(find(input(i).u(:,2) == 0.1), find(strcmp(tbl.phase, phase)));
                        elseif strcmp(opts{j}, 'medium')
                            idx = intersect(find(input(i).u(:,2) == 0.3), find(strcmp(tbl.phase, phase)));
                        elseif strcmp(opts{j}, 'difficult')
                            idx = intersect(find(input(i).u(:,2) == 0.9), find(strcmp(tbl.phase, phase)));
                        end
                    else
                        idx = intersect(find(strcmp(tbl.(field{1}), opts{j})), find(strcmp(tbl.phase, phase)));                    
                    end
                    subs(count) = i;
                    what(count) = 0; % logY
                    cond{count} = opts{j};
                    data(count) = mean(res.est(m,i).y(idx), 'omitnan');
                    count = count + 1;
                    subs(count) = i;
                    what(count) = 1; % logY
                    cond{count} = opts{j};
                    data(count) = mean(res.est(m,i).optim.yhat(idx), 'omitnan');
                    count = count + 1;
                end
            end
            temp.subID = subs;
            temp.cond  = categorical(cond);
            temp.real  = categorical(what,[0 1],{'logY','logYhat'});
            temp.logRT = data;
            out.(field{1}) = temp;
            fout{p} = out;
        end
    end
    % line graphs
    linesize = 3;
    figure
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.5 1]);
    count = 1;
    % comparison of the difficulty levels
    x = 1:3;
    for p = 1:3
        subplot(2,3,count)
        out = fout{p};
        count = count + 1;
        y    = [mean(out.difficulty(out.difficulty.cond == "easy" & out.difficulty.real == "logY", :).logRT),...
            mean(out.difficulty(out.difficulty.cond == "medium" & out.difficulty.real == "logY", :).logRT),...
            mean(out.difficulty(out.difficulty.cond == "difficult" & out.difficulty.real == "logY", :).logRT)];
        yerr = [std(out.difficulty(out.difficulty.cond == "easy" & out.difficulty.real == "logY", :).logRT),...
            std(out.difficulty(out.difficulty.cond == "medium" & out.difficulty.real == "logY", :).logRT),...
            std(out.difficulty(out.difficulty.cond == "difficult" & out.difficulty.real == "logY", :).logRT)];
        yerr = yerr/sqrt(nSubs);
        errorbar(x,y,yerr, 'LineWidth', linesize)
        hold on
        yhat    = [mean(out.difficulty(out.difficulty.cond == "easy" & out.difficulty.real == "logYhat", :).logRT),...
            mean(out.difficulty(out.difficulty.cond == "medium" & out.difficulty.real == "logYhat", :).logRT),...
            mean(out.difficulty(out.difficulty.cond == "difficult" & out.difficulty.real == "logYhat", :).logRT)];
        yhaterr = [std(out.difficulty(out.difficulty.cond == "easy" & out.difficulty.real == "logYhat", :).logRT),...
            std(out.difficulty(out.difficulty.cond == "medium" & out.difficulty.real == "logYhat", :).logRT),...
            std(out.difficulty(out.difficulty.cond == "difficult" & out.difficulty.real == "logYhat", :).logRT)];
        yhaterr = yhaterr/sqrt(nSubs);
        errorbar(x,yhat,yhaterr, 'LineWidth', linesize)
        xlim([0, 4])
        xticks([1, 2, 3])
        xticklabels({'easy', 'medium', 'difficult'})
        ylabel('logRT [ms]', 'FontSize', 14)
        xlabel(phases{p}, 'FontSize', 14)
        if p == 1
            legend('y', 'yhat')
        end
    end
    % comparison of the expected levels
    x = 1:2;
    for p = 1:3
        subplot(2,3,count)
        out = fout{p};
        count = count + 1;
        y    = [mean(out.expected(out.expected.cond == "expected" & out.expected.real == "logY", :).logRT),...
            mean(out.expected(out.expected.cond == "unexpected" & out.expected.real == "logY", :).logRT)];
        yerr = [std(out.expected(out.expected.cond == "expected" & out.expected.real == "logY", :).logRT),...
            std(out.expected(out.expected.cond == "unexpected" & out.expected.real == "logY", :).logRT)];
        yerr = yerr/sqrt(nSubs);
        errorbar(x,y,yerr, 'LineWidth', linesize)
        hold on
        yhat    = [mean(out.expected(out.expected.cond == "expected" & out.expected.real == "logYhat", :).logRT),...
            mean(out.expected(out.expected.cond == "unexpected" & out.expected.real == "logYhat", :).logRT)];
        yhaterr = [std(out.expected(out.expected.cond == "expected" & out.expected.real == "logYhat", :).logRT),...
            std(out.expected(out.expected.cond == "unexpected" & out.expected.real == "logYhat", :).logRT)];
        yhaterr = yhaterr/sqrt(nSubs);
        errorbar(x,yhat,yhaterr, 'LineWidth', linesize)
        xlim([0, 3])
        xticks([1, 2])
        xticklabels({'expected', 'unexpected'})
        xlabel(phases{p}, 'FontSize', 14)
        ylabel('logRT [ms]', 'FontSize', 14)
    end
    % save the plots
    print(strcat(figdir, filesep, 'aggLogRT-phases_model_', modSpace(m).name), '-dpng');
    print(strcat(figdir, filesep, 'aggLogRT-phases_model_', modSpace(m).name), '-dsvg');
    close;        
end

end

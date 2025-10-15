close all
clear
clc
dir_PAL  = pwd;                          
dir_HGF  = [dir_PAL filesep 'EMBA_HGF']; 
dir_data = [dir_PAL filesep 'data'];     
saveDir = [pwd filesep 'HGF_results_initialPriors'];

law = load('/home/emba/Documents/MATLAB/Propranolol-master/data/pupil_and_hgf_est/2/HGF_est.mat', 'est');

if ~exist(fullfile(saveDir), 'dir')
   mkdir(fullfile(saveDir))
end

% Pilot data set 
pilot = load("PAL_Prop_pilot.mat");
sdata = pilot.data(1:37);
load('PAL_data.mat', 'data');
sdata = data;
res   = struct();

% options.opt_config = 'tapas_quasinewton_optim_config';
% options.opt_config.maxStep = 2;

for n=1:length(sdata)

    sub = sdata(n);
    options = setupOpt_and_randSeed();

    % Invert the model
    est = tapas_fitModel(sub.y,... % responses
                sub.u,... % input sequence
                'emba_ehgf_binary_pu_tbt_config',...   %Prc fitting model
                'emba_logrt_linear_binary_C_config',... %Obs fitting model
                options.opt_config);                   %opt algo

    % add to res structure
    res.est(n) = est;
    
    % Save model fit as struct
    save_path = fullfile(saveDir, ['sub' num2str(n)], "est_mod");
    if ~exist(fullfile(saveDir, ['sub' num2str(n)]), 'dir')
       mkdir(fullfile(saveDir, ['sub' num2str(n)]))
    end
    save(save_path, '-struct', 'est');
end

save_path = fullfile(saveDir, 'full_results');
save(save_path, 'res');

%% Plot

nTrials = numel(res.est(1,1).y);
nSubs   = size(res.est, 2);

figdir = fullfile(saveDir, 'figures'); 
if ~exist(figdir, 'dir')
    mkdir(figdir)
end

% load stimulus order
tbl = readtable("PAL_scheme.csv");
coi = {'difficulty', 'expected', 'phase'}; % conditions of interest

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
            % find the corrisponding indices with difficulty being in
            % the u itself > different in Lawson and our study
            if strcmp(field{1}, 'difficulty')
                if strcmp(opts{j}, 'easy')
                    idx = find(res.est(i).u(:,2) == 0.1);
                elseif strcmp(opts{j}, 'medium')
                    idx = find(res.est(i).u(:,2) == 0.3);
                elseif strcmp(opts{j}, 'difficult')
                    idx = find(res.est(i).u(:,2) == 0.9);
                end
            else
                idx = find(strcmp(tbl.(field{1}), opts{j}));                    
            end
            subs(count) = i;
            what(count) = 0; % logY
            cond{count} = opts{j};
            data(count) = mean(res.est(i).y(idx), 'omitnan');
            count = count + 1;
            subs(count) = i;
            what(count) = 1; % logY
            cond{count} = opts{j};
            data(count) = mean(res.est(i).optim.yhat(idx), 'omitnan');
            count = count + 1;
        end
    end
    temp.subID = subs;
    temp.cond  = categorical(cond);
    temp.real  = categorical(what,[0 1],{'logY','logYhat'});
    temp.logRT = data;
    out.(field{1}) = temp;
end
% line graphs
figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.5 1]);
% comparison of the difficulty levels
subplot(3,1,1)
x = 1:3;
y    = [mean(out.difficulty(out.difficulty.cond == "easy" & out.difficulty.real == "logY", :).logRT),...
    mean(out.difficulty(out.difficulty.cond == "medium" & out.difficulty.real == "logY", :).logRT),...
    mean(out.difficulty(out.difficulty.cond == "difficult" & out.difficulty.real == "logY", :).logRT)];
yerr = [std(out.difficulty(out.difficulty.cond == "easy" & out.difficulty.real == "logY", :).logRT),...
    std(out.difficulty(out.difficulty.cond == "medium" & out.difficulty.real == "logY", :).logRT),...
    std(out.difficulty(out.difficulty.cond == "difficult" & out.difficulty.real == "logY", :).logRT)];
yerr = yerr/sqrt(nSubs);
errorbar(x,y,yerr)
hold on
yhat    = [mean(out.difficulty(out.difficulty.cond == "easy" & out.difficulty.real == "logYhat", :).logRT),...
    mean(out.difficulty(out.difficulty.cond == "medium" & out.difficulty.real == "logYhat", :).logRT),...
    mean(out.difficulty(out.difficulty.cond == "difficult" & out.difficulty.real == "logYhat", :).logRT)];
yhaterr = [std(out.difficulty(out.difficulty.cond == "easy" & out.difficulty.real == "logYhat", :).logRT),...
    std(out.difficulty(out.difficulty.cond == "medium" & out.difficulty.real == "logYhat", :).logRT),...
    std(out.difficulty(out.difficulty.cond == "difficult" & out.difficulty.real == "logYhat", :).logRT)];
yhaterr = yhaterr/sqrt(nSubs);
errorbar(x,yhat,yhaterr)
xlim([0, 4])
xticks([1, 2, 3])
xticklabels({'easy', 'medium', 'difficult'})
ylabel('logRT [ms]', 'FontSize', 14)
title('Comparison of subject-specific means of y and yhat', 'FontSize', 20)
legend('y', 'yhat')
% comparison of the expected levels
subplot(3,1,2)
x = 1:2;
y    = [mean(out.expected(out.expected.cond == "expected" & out.expected.real == "logY", :).logRT),...
    mean(out.expected(out.expected.cond == "unexpected" & out.expected.real == "logY", :).logRT)];
yerr = [std(out.expected(out.expected.cond == "expected" & out.expected.real == "logY", :).logRT),...
    std(out.expected(out.expected.cond == "unexpected" & out.expected.real == "logY", :).logRT)];
yerr = yerr/sqrt(nSubs);
errorbar(x,y,yerr)
hold on
yhat    = [mean(out.expected(out.expected.cond == "expected" & out.expected.real == "logYhat", :).logRT),...
    mean(out.expected(out.expected.cond == "unexpected" & out.expected.real == "logYhat", :).logRT)];
yhaterr = [std(out.expected(out.expected.cond == "expected" & out.expected.real == "logYhat", :).logRT),...
    std(out.expected(out.expected.cond == "unexpected" & out.expected.real == "logYhat", :).logRT)];
yhaterr = yhaterr/sqrt(nSubs);
errorbar(x,yhat,yhaterr)
xlim([0, 3])
xticks([1, 2])
xticklabels({'expected', 'unexpected'})
ylabel('logRT [ms]', 'FontSize', 14)
% comparison of the phases
subplot(3,1,3)
x = 1:3;
y    = [mean(out.phase(out.phase.cond == "prevolatile" & out.phase.real == "logY", :).logRT),...
    mean(out.phase(out.phase.cond == "volatile" & out.phase.real == "logY", :).logRT),...
    mean(out.phase(out.phase.cond == "postvolatile" & out.phase.real == "logY", :).logRT)];
yerr = [std(out.phase(out.phase.cond == "prevolatile" & out.phase.real == "logY", :).logRT),...
    std(out.phase(out.phase.cond == "volatile" & out.phase.real == "logY", :).logRT),...
    std(out.phase(out.phase.cond == "postvolatile" & out.phase.real == "logY", :).logRT)];
yerr = yerr/sqrt(nSubs);
errorbar(x,y,yerr)
hold on
yhat    = [mean(out.phase(out.phase.cond == "prevolatile" & out.phase.real == "logYhat", :).logRT),...
    mean(out.phase(out.phase.cond == "volatile" & out.phase.real == "logYhat", :).logRT),...
    mean(out.phase(out.phase.cond == "postvolatile" & out.phase.real == "logYhat", :).logRT)];
yhaterr = [std(out.phase(out.phase.cond == "prevolatile" & out.phase.real == "logYhat", :).logRT),...
    std(out.phase(out.phase.cond == "volatile" & out.phase.real == "logYhat", :).logRT),...
    std(out.phase(out.phase.cond == "postvolatile" & out.phase.real == "logYhat", :).logRT)];
yhaterr = yhaterr/sqrt(nSubs);
errorbar(x,yhat,yhaterr)
xlim([0, 4])
xticks([1, 2, 3])
xticklabels({'prevolatile', 'volatile', 'postvolatile'})
ylabel('logRT [ms]', 'FontSize', 14)
% save the plots
print(strcat(figdir, filesep, 'aggLogRT_model'), '-dpng');
% close;        

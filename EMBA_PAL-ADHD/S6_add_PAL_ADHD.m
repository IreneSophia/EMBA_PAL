%% run a simpler model, specifically the Rescorla-Wagner learning model

close all
clear
clc

%% Specify the necessary parameters

% Current directory to be able to return to it after toolbox installation
dir_PAL  = pwd;                          

% Path to the EMBA_HGF toolbox
dir_HGF  = [fileparts(dir_PAL) filesep 'EMBA_HGF']; 

% Path to all the data files
dir_data = [fileparts(dir_PAL) filesep 'data'];     

% The directory, where the results will be saved.
saveDir = [pwd filesep 'HGF_results'];
addDir  = [saveDir filesep 'addRW'];

% Create the directory for the revision model
if ~exist(addDir, 'dir')
    mkdir(addDir)
end

% Number of simulations (number of simulated subjects)
nSim = 100; 

% Number of local cores for parallelisation
local_cores = 32;

% Switch for parallelisation: 1 for parallel, 0 for the serial version
if ismac()
    parallel = 0;
else
    parallel = 1;
end

%% Setup the toolboxes

addpath([dir_HGF filesep '/tapas-master']);
tapas_init('hgf');

if ~ismac
    cd([dir_HGF filesep 'VBA-toolbox-master']); 
    VBA_setup();
    cd(dir_PAL)
end

%% Add the relevant subfolders to the path

addpath(dir_HGF);
addpath(dir_data);
addpath([dir_HGF filesep 'models']);
addpath([dir_HGF filesep 'workflow']);

%% Pipeline and model space setup

% Get the observation and perception model names
obsNamesDisp = "SUR";
prcNamesDisp = "RW";
prcNamesList = "tapas_rw_binary";
obsNamesList = "emba_logrt_linear_binary_SUR";

% Create a structure with the model space
modSpaceNew = setupModelSpace(obsNamesDisp, prcNamesDisp, obsNamesList, prcNamesList, true);

% Add the HGF models to the model space
[obsNamesDisp, prcNamesDisp, obsNamesList, prcNamesList] = modelNames_PAL_ADHD();
modSpaceHGF = setupModelSpace(obsNamesDisp, prcNamesDisp, obsNamesList, prcNamesList, true);

% Load the participant data
load('PAL-ADHD_data.mat', 'data'); 

%% Model inversion of pilot data

% Pilot data set 
pilot = load("PAL_Prop_pilot.mat");

% Number of participants
nPilots = numel(pilot.data);

% fit the model - takes a bit
fitPilotData(modSpaceNew, nPilots, pilot, parallel, local_cores, addDir)

% Accumulate the results. Necessary for some of the plotting. 
accumulateForPlotting('pilots', modSpaceNew, addDir)

%% Quick posterior predictive checks for the pilots

% Plot predicted differences between difficulty & expectedness [!HARDCODED!]
plotAggLogRTsModel(modSpaceNew, 'pilots', addDir, false, true, pilot.data) % [!ADJUSTED]

%% Empirical priors computation

% Compute and save the priors
computeEmpiricalPriors(modSpaceNew, nPilots, 'pilots', addDir);

% Plot empirical priors [!HARDCODED!]
plotEmpiricalPriors(modSpaceNew, 'pilots', addDir, []); 

% Update the model space to use empirical priors 
modSpaceRobMean = updateModSpace(modSpaceNew, addDir);

% Second version of the RW model with default priors
modSpaceRobMean = [modSpaceRobMean modSpaceNew];
modSpaceRobMean(end).name = modSpaceRobMean(end).name + "default";
modSpaceRobMean(end).prc_config.logitalsa   = 1;
modSpaceRobMean(end).prc_config.priorsas(2) = 1;

% Combine with updated HGF space
modSpaceAll = [updateModSpace(modSpaceHGF, saveDir) modSpaceRobMean];

%% Simulate the data 
% Here, we assume that the inputs u are the same for everyone
% We can also decide how verbose the output is supposed to be and the
% number of maximum tries in case of instability
verbose = false;
maxRep  = 10;
simModelsAppend(modSpaceRobMean, nSim, pilot.data(end).u, saveDir, addDir, verbose, maxRep);

%% Fit all models to simulated data
% This now also needs to rerun the other models with the current one 
% to ensure that we can do the Model Identifiability across the models. 

fitSimDataAppend(modSpaceAll, nSim, parallel, local_cores, saveDir, addDir);

%% Recovery analysis (parameter recovery + model identifyability) 

checkParameterRecovery(modSpaceAll, nSim, addDir);
plotParameterRecovery(modSpaceAll, addDir); % [!HARDCODED!]

% Here, we see that the empirical priors are too narrow for the RW model.
% We adjusted the width of the prior above to account for that and now
% disregard the model with the narrow, empirical priors. 
modSpaceSel = modSpaceAll([1:4,end]);

% Save all the model spaces for the future
save(fullfile(addDir, 'modSpace.mat'), 'modSpaceSel', 'modSpaceAll',...
    'modSpaceRobMean', 'modSpaceHGF', 'modSpaceNew')

% now we check the model identifiability
checkModelIdentifiability(nSim, modSpaceSel, addDir);
plotModelIdentifiability(modSpaceSel, addDir);

%% Fit the models on the data

% Fit the models to the participant data - default RW-SUR model
fitData(modSpaceNew, data, parallel, local_cores, addDir)

% Accumulate the results. Necessary for some of the plotting. 
accumulateForPlotting('main', modSpaceNew, addDir)

%% Compare all models

% load HGF models
hgf = load([saveDir filesep 'main' filesep 'full_results.mat']);
new = load([addDir  filesep 'main' filesep 'full_results.mat']);

% combine
est = [hgf.res.est; new.res.est];

% extract the LME
LME = nan(size(est));
txt = cell(size(est,1),1);
for m = 1:size(est,1)
    txt{m} = modSpaceSel(m).name;
    for i = 1:size(est,2)
        LME(m,i) = est(m,i).optim.LME;
    end
end

% compare models over all participants
options.modelNames = txt;
options.verbose    = false;
[posterior, out]   = VBA_groupBMC(LME, options);

% save the outcome and plots
print([addDir filesep 'model_comparison'], '-dpng');
print([addDir filesep 'model_comparison'], '-dsvg');
close;
save([addDir filesep 'model_comparison.mat'], 'posterior', 'out')
writetable(table(out.Ef, out.ep.', 'VariableNames', {'Ef', 'ep'}), [addDir filesep 'model_comparison_out.csv']);
winModel = find(out.Ef == max(out.Ef)); % winModel = 2;
fprintf('\n\nThe winning model is %s!\n\n', out.options.modelNames{winModel})

%% Optimal parameters for the winning model

bopars = tapas_fitModel([],...
                        data(1).u,...
                        modSpaceSel(winModel).prc_config,...
                        'tapas_bayes_optimal_binary_config');

% extract values for plotting
tbl = table();
tbl.trl    = [1:size(bopars.traj.wt,1)].';
tbl.alpha1_opt = bopars.traj.wt(:,1);
tbl.alpha2_opt = bopars.traj.wt(:,2);
tbl.alpha3_opt = bopars.traj.wt(:,3);
writetable(tbl, fullfile(addDir, 'BOpars.csv'))


%% Perform posterior predictive checks for RW model

% Plot simulated data in comparison to real subject data for RW model
plotLogRTsModel(modSpaceNew, 'main', addDir, 0, 1);

% Create simulated data for better plotting in R
simPostRT(modSpaceNew, 'main', addDir, 1, nSim);

% Plot predicted differences between difficulty & expectedness [!HARDCODED!]
plotAggLogRTsModel(modSpaceNew, 'main', addDir, false, true, data)
plotAggLogRTsModelPhases(modSpaceHGF, 'main', saveDir, data)

% Plot average regressors of winning model [!HARDCODED]
plotAverageRegressors(modSpaceSel(winModel), 'main', data, saveDir);

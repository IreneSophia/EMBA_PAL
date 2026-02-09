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
revDir  = [saveDir filesep 'rev'];

% Create the directory for the revision model
if ~exist(revDir, 'dir')
    mkdir(revDir)
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
[obsNamesDisp, prcNamesDisp, obsNamesList, prcNamesList] = modelNames_PAL_ASD();
modSpaceHGF = setupModelSpace(obsNamesDisp, prcNamesDisp, obsNamesList, prcNamesList, false);

% Load the participant data
load('PAL-ASD_data.mat', 'data'); 

%% Model inversion of pilot data

% Pilot data set 
pilot = load("PAL_Prop_pilot.mat");

% Number of participants
nPilots = numel(pilot.data);

% fit the model - takes a bit
fitPilotData(modSpaceNew, nPilots, pilot, parallel, local_cores, revDir)

% Accumulate the results. Necessary for some of the plotting. 
accumulateForPlotting('pilots', modSpaceNew, revDir)

%% Quick posterior predictive checks for the pilots

% Plot predicted differences between difficulty & expectedness [!HARDCODED!]
plotAggLogRTsModel(modSpaceNew, 'pilots', revDir, false, true, pilot.data(end).u) % [!ADJUSTED]

%% Empirical priors computation

% Compute and save the priors
computeEmpiricalPriors(modSpaceNew, nPilots, 'pilots', revDir);

% Plot empirical priors [!HARDCODED!]
plotEmpiricalPriors(modSpaceNew, 'pilots', revDir, []); 

% Update the model space to use empirical priors 
modSpaceRobMean = updateModSpace(modSpaceNew, revDir);

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
simModelsAppend(modSpaceRobMean, nSim, pilot.data(end).u, saveDir, revDir, verbose, maxRep);

%% Fit all models to simulated data
% This now also needs to rerun the other models with the current one 
% to ensure that we can do the Model Identifiability across the models. 

fitSimDataAppend(modSpaceAll, nSim, parallel, local_cores, saveDir, revDir);

%% Recovery analysis (parameter recovery + model identifyability) 

checkParameterRecovery(modSpaceAll, nSim, revDir);
plotParameterRecovery(modSpaceAll, revDir); % [!HARDCODED!]

% Here, we see that the empirical priors are too narrow for the RW model.
% We adjusted the width of the prior above to account for that and now
% disregard the model with the narrow, empirical priors. 
modSpaceSel = modSpaceAll([1,2,end]);
modSpaceRobMean = modSpaceNew;

% now we check the model identifiability
checkModelIdentifiability(nSim, modSpaceSel, revDir); % *** VBA_random: alpha must be a vector of positive values.
plotModelIdentifiability(modSpaceAll, revDir);

%% Fit the models on the data

% Fit the models to the participant data
fitData(modSpaceRobMean, data, parallel, local_cores, revDir)

% Accumulate the results. Necessary for some of the plotting. 
accumulateForPlotting('main', modSpaceRobMean, revDir)

%% Compare all models

% load HGF models
hgf = load([saveDir filesep 'main' filesep 'full_results.mat']);
new = load([revDir  filesep 'main' filesep 'full_results.mat']);

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
print([revDir filesep 'model_comparison'], '-dpng');
print([revDir filesep 'model_comparison'], '-dsvg');
close;
save([revDir filesep 'model_comparison.mat'], 'posterior', 'out')
writetable(table(out.Ef, out.ep.', 'VariableNames', {'Ef', 'ep'}), [revDir filesep 'model_comparison_out.csv']);
winModel = find(out.Ef == max(out.Ef)); % winModel = 2;
fprintf('\n\nThe winning model is %s!\n\n', out.options.modelNames{winModel})

%% Perform posterior predictive checks

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
parallel = 0;

%% Setup the toolboxes

addpath([dir_HGF filesep '/tapas-master']);
tapas_init('hgf');

% cd([dir_HGF filesep 'VBA-toolbox-master']); 
% VBA_setup();
% cd(dir_PAL)

%% Add the relevant subfolders to the path

addpath(dir_HGF);
addpath(dir_data);
addpath([dir_HGF filesep 'models']);
addpath([dir_HGF filesep 'workflow']);

%% Pipeline and model space setup

% Get the observation and perception model names
obsNamesDisp = "GAU";
prcNamesDisp = "RW";
prcNamesList = "tapas_rw_binary";
obsNamesList = "tapas_gaussian_obs";

% Create a structure with the model space
modSpaceNew = setupModelSpace(obsNamesDisp, prcNamesDisp, obsNamesList, prcNamesList, true);

% Add the HGF models to the model space
[obsNamesDisp, prcNamesDisp, obsNamesList, prcNamesList] = modelNames_PAL_ASD();
modSpaceHGF = setupModelSpace(obsNamesDisp, prcNamesDisp, obsNamesList, prcNamesList, false);
modSpaceAll = [updateModSpace(modSpaceHGF, saveDir) modSpaceNew];

% Load the participant data
load('PAL-ASD_data.mat', 'data'); 

%% Simulate the data 
% Here, we assume that the inputs u are the same for everyone
% We can also decide how verbose the output is supposed to be and the
% number of maximum tries in case of instability
verbose = false;
maxRep  = 10;
simModelsAppend(modSpaceNew, nSim, data(end).u, saveDir, revDir, verbose, maxRep);

%% Fit all models to simulated data
% This now also needs to rerun the other models with the current one 
% to ensure that we can do the Model Identifiability across the models. 

fitSimDataAppend(modSpaceAll, nSim, parallel, local_cores, saveDir, revDir);

%% Recovery analysis (parameter recovery + model identifyability) [!MISSING]

checkParameterRecovery(modSpaceRobMeanAll, nSim, revDir);
plotParameterRecovery(modSpaceRobMeanAll, revDir); % [!HARDCODED!]

checkModelIdentifiability(nSim, modSpaceRobMeanAll, revDir); %
plotModelIdentifiability(modSpaceRobMeanAll, revDir);

%% Fit the models on the data

% Fit the models to the participant data
fitData(modSpaceNew, data, parallel, local_cores, revDir)

% Accumulate the results. Necessary for some of the plotting. 
accumulateForPlotting('main', modSpaceNew, revDir)


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
    txt{m} = modSpaceAll(m).name;
    for i = 1:size(est,2)
        LME(m,i) = est(m,i).optim.LME;
    end
end

% compare models 
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
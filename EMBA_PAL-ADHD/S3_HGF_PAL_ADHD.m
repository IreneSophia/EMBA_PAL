% Main running file of the EMBA_HGF toolbox
% Before running the script
%     - navigate to the directory where this script is located
%     - specify the models in the corresponding getModelNames file in this
%     folder (example is in the toolbox)
% 
% CAUTION : some of the plotting is hardcoded for the EMBA models

%-----------------------------------------------------------------------------
% 
% Copyright (C) 2024 Anna Yurova, Irene Sophia Plank, LMU University Hospital; 
%               based on the hessetal_spirl_analysis toolbox by Alex Hess (2024), TNU, ETHZ
%
% Please cite the papers from CITATION.txt when using the EMBA_HGF Toolbox.
%
% This file is part of the EMBA_HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. For further details, see the file
% LICENSE or <http://www.gnu.org/licenses/>.
%-----------------------------------------------------------------------------

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

% Number of simulations (number of simulated subjects)
nSim = 100; 

% Number of local cores for parallelisation
local_cores = 32;

% Switch for parallelisation: 1 for parallel, 0 for the serial version
parallel = 0;

%% Setup the toolboxes

addpath([dir_HGF filesep '/tapas-master']);
tapas_init('hgf');

cd([dir_HGF filesep 'VBA-toolbox-master']); 
VBA_setup();
cd(dir_PAL)

%% Add the relevant subfolders to the path

addpath(dir_HGF);
addpath(dir_data);
addpath([dir_HGF filesep 'models']);
addpath([dir_HGF filesep 'workflow']);

%% Pipeline setup

% Get the observation and perception model names
[obsNamesDisp, prcNamesDisp, obsNamesList, prcNamesList] = modelNames_PAL_ADHD();

% Pilot data set 
pilot = load("PAL_Prop_pilot.mat");

% Number of participants
nPilots = numel(pilot.data);

%% Model space setup

% Create a structure with the model space
modSpace = setupModelSpace(obsNamesDisp, prcNamesDisp, obsNamesList, prcNamesList, true);

% Plot initial priors for the models [!HARDCODED!]
plotInitialPriorsPredDens(modSpace, pilot.data(end).u, saveDir)

%% Model inversion
% Fit the all models for each pilot subject
% The results are saved in .mat files. This line can be commented out once
% the fitting is completed. 
% Do not forget to uncomment if any setup has changed!

fitPilotData(modSpace, nPilots, pilot, parallel, local_cores, saveDir)

% Accumulate the results. Necessary for some of the plotting. 
accumulateForPlotting('pilots', modSpace, saveDir)

%% Quick posterior predictive checks for the pilots

% Plot predicted differences between difficulty & expectedness [!HARDCODED!]
plotAggLogRTsModel(modSpace, 'pilots', saveDir, false, true)

% Plot regressors [!HARDCODED!]
plotAverageRegressors(modSpace, 'pilots', pilot.data, saveDir);

%% Empirical priors computation

% Compute and save the priors
computeEmpiricalPriors(modSpace, nPilots, 'pilots', saveDir);

% Plot empirical priors [!HARDCODED!]
plotEmpiricalPriors(modSpace, 'pilots', saveDir, []);

% Update the model space to use empirical priors
modSpaceRobMean = updateModSpace(modSpace, saveDir);

%% Simulate the data with the empirical priors
% Here, we assume that the inputs u are the same for everyone
% We can also decide how verbose the output is supposed to be and the
% number of maximum tries in case of instability
verbose = false; 
maxRep  = 10;
simModels(modSpaceRobMean, nSim, pilot.data(end).u, saveDir, verbose, maxRep);
plotSimData(modSpaceRobMean, nSim, saveDir); %[!HARDCODED!]

%% Fit all models to simulated data
% The results are saved in .mat files. This line can be commented out once
% the fitting is completed. 
% Do not forget to uncomment if any setup has changed!

fitSimData(modSpaceRobMean, nSim, parallel, local_cores, saveDir); % [!CURRENTLY HERE!]

%% Recovery analysis (parameter recovery + model identifyability)

checkParameterRecovery(modSpaceRobMean, nSim, saveDir);
plotParameterRecovery(modSpaceRobMean, saveDir); % [!HARDCODED!]

checkModelIdentifiability(nSim, modSpaceRobMean, saveDir); %
plotModelIdentifiability(modSpaceRobMean, saveDir);

%% Plot the simulations

% Plot RTs of the simulations [!HARDCODED!]
plotSimRec(modSpace, nSim, saveDir); 
plotLogRTsModel(modSpace, 'sim', saveDir, 8, []); % [!UNTIL HERE]

%% Fit the models on the ACTUAL DATA

% Load the participant data
load('PAL-ADHD_data.mat', 'data'); 

% Fit the models to the participant data
fitData(modSpaceRobMean, data, parallel, local_cores, saveDir)

% Accumulate the results. Necessary for some of the plotting. 
accumulateForPlotting('main', modSpaceRobMean, saveDir)

%% Compare models 

% Similar to model identifiability > saves the figures
out = compareModels(modSpaceRobMean, saveDir); %
save([saveDir filesep 'main' filesep 'model_comparison_out'], 'out');
writetable(table(out.Ef, out.ep.', 'VariableNames', {'Ef', 'ep'}), 'model_comparison_out.csv');
winModel = find(out.Ef == max(out.Ef)); %winModel = 2;
fprintf('\n\nThe winning model is %s!\n\n', out.options.modelNames{winModel})

%% Perform posterior predictive checks

% Plot trajectoris for single participants for winning model [!HARDCODED!]
plotTrajForMain(modSpaceRobMean, winModel, parallel, local_cores, saveDir);

% Plot simulated data in comparison to real subject data for winning model
plotLogRTsModel(modSpaceRobMean, 'main', saveDir, 0, winModel);

% Plot many simulated datasets in comparison to real subject data for winning model
% focusing on the best and the worst fit participants 
plotMultiLogRTsModel(modSpaceRobMean, 'main', saveDir, ...
    [67, 64, 38, 71, 19, 82, 36, 73], winModel, nSim, true);

% Create simulated data for better plotting in R
simPostRT(modSpaceRobMean, 'main', saveDir, winModel, nSim);

% Plot predicted differences between difficulty & expectedness [!HARDCODED!]
plotAggLogRTsModel(modSpaceRobMean(winModel), 'main', saveDir, true, true)
plotVarLogRTsModel(modSpaceRobMean(winModel), 'main', saveDir, true, true)

% Plot empirical priors with parameters extracted from data [!HARDCODED!]
computeEmpiricalPriors(modSpaceRobMean(winModel), length(data), 'main', saveDir);
plotEmpiricalPriors(modSpaceRobMean(winModel), 'main', saveDir, []);

%% Extract model parameters for further analysis

% extrac the "raw" values from the winning model
extractHGFcsv(modSpaceRobMean, saveDir, winModel)

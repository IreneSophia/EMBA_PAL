% Main running file of the EMBA_HGF toolbox
% Before running the script
%     - navigate to the directory where the "EMBA_PAL_main" script is located
%     - specify the models in the getModelNames file
% 
% CAUTION : some of the plotting is hardcoded to the EMBA models

%-----------------------------------------------------------------------------
% 
% Copyright (C) 2024 Anna Yurova, Irene Sophia Plank, LMU University Hospital; 
%               based on the hessetal_spirl_analysis toolbox by Alex Hess (2024), TNU, ETHZ
%
% Please cite <INSERT REFERENCE>, Hess et. al. (2024) from
% CITATION.txt when using the EMBA_HGF Toolbox.
%
% This file is part of the EMBA_HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. For further details, see the file
% LICENSE or <http://www.gnu.org/licenses/>.
%-----------------------------------------------------------------------------

close all
clear
clc

dir_PAL = pwd;
dir_HGF = [dir_PAL filesep 'EMBA_HGF'];

%% Setup the toolboxes

addpath([dir_HGF filesep '/tapas-master']);
tapas_init('hgf');

cd([dir_HGF filesep 'VBA-toolbox-master']); 
VBA_setup();
cd(dir_PAL)

addpath([dir_HGF filesep 'RainCloudPlots-master/tutorial_matlab/']);

%% Add the relevant subfolders to the path

addpath(dir_HGF);
addpath([dir_HGF filesep 'data']);
addpath([dir_HGF filesep 'models']);
addpath([dir_HGF filesep 'workflow']);

%% Specify the necessary parameters

% Pilot data set 
pilot = load("PAL_Prop_pilot.mat");

% The directory, where the results will be saved.
saveDir = [pwd filesep 'HGF_results'];

% Number of simulations (number of simulated subjects)
nSim = 100; 

%% Parallel toolbox setup

% Number of local cores
local_cores = 32;

% Switch for parallelization: 1 for parallel version, 0 for the serial
% version
parallel = 1;

%% Pipeline setup

% Get the observation and perception model names
[obsNamesDisp, prcNamesDisp, obsNamesList, prcNamesList] = modelNames();

% Number of participants
nPilots = numel(pilot.data);

%% Model space setup
modSpace = setupModelSpace(obsNamesDisp, prcNamesDisp, obsNamesList, prcNamesList);

% Number of models
nModels = size(modSpace, 2);

% Plot initial priors for the models [!HARDCODED!]
plotInitialPriorsPredDens(modSpace, pilot.data(end).u, saveDir)

%% Model inversion
% Fit the all models for each pilot subject
% The results are saved in .mat files. This line can be commented out once
% the fitting is completed. 
% Do not forget to uncomment if any setup has changed!

fitPilotData(modSpace, nPilots, pilot, parallel, local_cores, saveDir)

% Plot regressors [!HARDCODED!] [! NEEDS ADJUSTING !]
plotAverageRegressors(modSpace, 'pilots', pilot.data, saveDir);

% Plot RTs with yhats
plotLogRTsModel(modSpace, 'pilots', saveDir); 

%% Empirical priors computation

% Compute and save the priors
computeEmpiricalPriors(modSpace, nPilots, saveDir);

% Plot empirical priors [!HARDCODED!]
plotEmpiricalPriors(modSpace, saveDir);

% Update the model space to use empirical priors
modSpaceRobMean = updateModSpace(modSpace, saveDir);

%% Simulate the data with the empirical priors
% Here, we assume that the inputs u are the same for everyone
% We can also decide how verbose the output is supposed to be and the
% number of maximum tries in case of instability
verbose = false;
maxRep  = 10;
simModels(modSpaceRobMean, nSim, pilot.data(end).u, saveDir, verbose, maxRep);
plotSimData(modSpaceRobMean, nSim, saveDir);

%% Fit all models to simulated data
% The results are saved in .mat files. This line can be commented out once
% the fitting is completed. 
% Do not forget to uncomment if any setup has changed!

fitSimData(modSpaceRobMean, nSim, parallel, local_cores,  saveDir);

%% Recovery analysis (parameter recovery + model identifyability)
checkParameterRecovery(modSpaceRobMean, nSim, saveDir);
plotParameterRecovery(modSpaceRobMean, saveDir); % [!parameter names are HARDCODED!]

checkModelIdentifiability(nSim, modSpaceRobMean, saveDir); %
plotModelIdentifiability(modSpaceRobMean, saveDir);

%% Plot the simulations

% Plot RTs of the simulations
plotSimRec(modSpace, nSim, saveDir); 

% Plot trajectories based on all of the simulations
plotTrajForAllSim(modSpace, nSim, saveDir);   % [! NEEDS ADJUSTING !]

%% Fit the preregistered model on the ACTUAL DATA

% Reduce model space to preregistered model
modSpaceWin = modSpaceRobMean(1);

% Load the participant data
load('PAL_data.mat', 'data');

% Fit the winning model to the participant data
fitData(modSpaceWin, data, parallel, local_cores, saveDir)

%% Perform posterior predictive checks


%% Extract the model parameters for further analysis

% Extract the relevant parameter and save them to a csv file
extractHGFcsv(modSpaceWin, saveDir)

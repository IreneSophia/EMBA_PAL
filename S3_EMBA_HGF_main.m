% Main running file of the EMBA_HGF toolbox
% Before running the script
%     - navigate to the directory where the "EMBA_HGF_main" script is located
%     - specify the models in the getModelNames file
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
dir_HGF  = [dir_PAL filesep 'EMBA_HGF']; 

% Path to all the data files
dir_data = [dir_PAL filesep 'data'];     

% The directory, where the results will be saved.
saveDir = [pwd filesep 'HGF_results'];

% Number of simulations (number of simulated subjects)
nSim = 100; 

% Number of local cores for parallelisation
local_cores = 32;

% Switch for parallelisation: 1 for parallel, 0 for the serial version
parallel = 1;

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
[obsNamesDisp, prcNamesDisp, obsNamesList, prcNamesList] = modelNames();

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

% We only keep models that actually capture two things: the difference
% between difficult, medium and easy trials as well as the difference
% between expected and unexpected trials. This applies to all models with
% the response model including the precision-weighted prediction error on
% level 2, coded as C. Therefore, we continue with those models. 
idx = contains(obsNamesDisp, "C");
obsNamesDisp = obsNamesDisp(idx);
obsNamesList = obsNamesList(idx);

% Create a structure with the updated model space
modSpace = setupModelSpace(obsNamesDisp, prcNamesDisp, obsNamesList, prcNamesList, true);

% Re-accumulate the results. Necessary for some of the plotting. 
accumulateForPlotting('pilots', modSpace, saveDir)

%% Continue with models passing posterior predictive check

% Plot regressors [!HARDCODED!]
plotAverageRegressors(modSpace, 'pilots', pilot.data, saveDir);

% Plot RTs with yhats [!HARDCODED!]
plotLogRTsModel(modSpace, 'pilots', saveDir, 0, []); 

%% Empirical priors computation

% Compute and save the priors
computeEmpiricalPriors(modSpace, nPilots, 'pilots', saveDir);

% Plot empirical priors [!HARDCODED!]
plotEmpiricalPriors(modSpace, 'pilots', saveDir);

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

% Some of the models produce suboptimal simulated data. We still continue
% with all of them to have a better idea of their performance. 

%% Fit all models to simulated data
% The results are saved in .mat files. This line can be commented out once
% the fitting is completed. 
% Do not forget to uncomment if any setup has changed!

fitSimData(modSpaceRobMean, nSim, parallel, local_cores,  saveDir);

%% Recovery analysis (parameter recovery + model identifyability)

checkParameterRecovery(modSpaceRobMean, nSim, saveDir);
plotParameterRecovery(modSpaceRobMean, saveDir); % [!HARDCODED!]

checkModelIdentifiability(nSim, modSpaceRobMean, saveDir); %
plotModelIdentifiability(modSpaceRobMean, saveDir);

%% Plot the simulations

% Plot RTs of the simulations [!HARDCODED!]
plotSimRec(modSpace, nSim, saveDir); 
plotLogRTsModel(modSpace, 'sim', saveDir, 8, []);

%% Fit the models on the ACTUAL DATA

% Load the participant data
load('PAL_data.mat', 'data');

% We fit all the models that passed the initial posterior predictive checks
% based on the pilot data. Additionally, we will add the model that was
% used by Lawson et al. (2017) with the priors that they used. 
[lawModel, modSpace] = addLawModelSpace(modSpace);
[~, modSpaceRobMean] = addLawModelSpace(modSpaceRobMean);

% Fit the models to the participant data
fitData(modSpaceRobMean, data, parallel, local_cores, saveDir)

% Accumulate the results. Necessary for some of the plotting. 
accumulateForPlotting('main', modSpace, saveDir)

%% Compare models 

% Similar to model identifiability > saves the figures
out = compareModels(modSpaceRobMean, saveDir); %
winModel = find(max(out.Ef));
fprintf('\n\nThe winning model is %s!\n\n', out.options.modelNames{winModel})

% This shows that the best model for the actual data is the eHGF-C which 
% also did not show any red flags for the simulated data. 

%% Perform posterior predictive checks

% Plot trajectoris for single participants for winning model [!HARDCODED!]
plotTrajForMain(modSpace, winModel, parallel, local_cores, saveDir);  

% Plot simulated data in comparison to real subject data for winning model
plotLogRTsModel(modSpace, 'main', saveDir, 0, winModel);

% Plot predicted differences between difficulty & expectedness [!HARDCODED!]
plotAggLogRTsModel(modSpace, 'main', saveDir, false, true)

% Plot empirical priors with parameters extracted from data [!HARDCODED!]
computeEmpiricalPriors(modSpace(1:end-1), length(data), 'main', saveDir);
plotEmpiricalPriors(modSpace, 'main', saveDir, winModel);

%% Extract the model parameters for further analysis

% We extract the parameters from the winning model to test our hypotheses
% and from the Lawson model for comparison
extractHGFcsv(modSpaceRobMean, saveDir, [winModel lawModel])

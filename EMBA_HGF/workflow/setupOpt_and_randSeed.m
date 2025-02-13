function options = setupOpt_and_randSeed()
% Specify the following options:
%   - optimization algorithm
%   - number of re-initializations (!!! affects performance !!!)
%   - random number generator options
%
% Note that this function resets the random number generator seed.
% Based on hessetal_spirl_analysis/utils/spirl_specs.m
%
%-----------------------------------------------------------------------------
%
% Copyright (C) 2024 Anna Yurova, Irene Sophia Plank, LMU University Hospital;
%               based on the hessetal_spirl_analysis toolbox by Alex Hess (2024), TNU, ETHZ
%
%
% This file is part of the EMBA_HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. For further details, see the file LICENSE or <http://www.gnu.org/licenses/>.
%-----------------------------------------------------------------------------

% Optimization algorithm
options.opt_config = eval('tapas_quasinewton_optim_config');

% Number of random re-initialisations (!!! affects performance !!!)
% Increase in case there are a lot of non-converging results
% In hessetal_spirl_analysis was set to 399
options.opt_config.nRandInit = 1;

% Set up random number generator
rng(123, 'twister');
options.rng.settings = rng;
options.rng.idx = 1; % Set counter for random number states

options.opt_config.seedRandInit = options.rng.settings.State(options.rng.idx, 1);

end

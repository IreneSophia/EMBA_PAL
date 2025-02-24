*Version 1.2*

EMBA analysis pipeline: Bayesian workflow for HGF modelling (EMBA_HGF)
========================================================================

------------
Contributors 
------------

|                                     |                                             |
| ----------------------------------- | ------------------------------------------- |
| Project lead, developer:            | Irene Sophia Plank                          |
| Computational scientist, developer: | Anna Yurova                                 |
| NEVIA lab leader:                   | Christine M. Falter-Wagner                  |
| Project members:                    | Alexandra Pior, Julia Nowak                 |
| Study:                              | EMBA (study acronym)                        |
| Date:                               | <INSERT DATE>                               |
| License:                            | GNU General Public License v3.0             |


-----------
DESCRIPTION
-----------

EMBA_HGF is part of the analysis pipeline for the EMBA study:

<INSERT REFERENCE>

The study was conducted at the NEVIA lab led by Christine M. Falter-Wagner at the LMU University Hospital in Munich. The project lead is Irene Sophia Plank.

The EMBA_HGF pipeline is based on the hessetal_spirl_analysis toolbox by Alex Hess (2024), TNU, ETHZ
which was introduced in 

Hess, A.J., Iglesias, S., Köchli, L., Marino, S., Müller-Schrader, M., 
Rigoux, L., Mathys, C., Harrison, O.K., Heinzle, J., Frässle, S., Stephan, K.E., 2024. 
Bayesian Workflow for Generative Modeling in Computational Psychiatry. 
*bioRxiv*, https://doi.org/10.1101/2024.02.19.581001

The first version of this pipeline was developed by Anna Yurova in a separate repository. The toolbox was moved to this repository after Irene Sophia Plank started adapting the pipeline. 

------------
INSTALLATION
------------

EMBA_HGF is written in MATLAB, make sure you have a valid installation of MATLAB.

Please download EMBA_HGF from the general PAL repository of the project:

https://github.com/IreneSophia/EMBA_PAL

Download and extract the following toolboxes into the main folder:
- [tapas](https://github.com/translationalneuromodeling/tapas/releases/tag/v6.0.2)
- [VBA-toolbox](https://mbb-team.github.io/VBA-toolbox/)
- [RainCloudPlots](https://github.com/RainCloudPlots/RainCloudPlots)

This should result in three new folders in the main directory: 
- tapas-master
- VBA-toolbox-master
- RainCloudPlots-master

--------
GET DATA
--------

Add your data to the "data" folder. 
The data should be stored as a ".mat" file and have the following structure:

- It should be a ( 1 x Number_of_Subjects ) struct
- Each struct entry should contain the following fields:
  - u: Number_of_Trials x 2 array of the inputs
  - y: Number_of_Trials     array of response times
  
This pipeline assumes that u is the same for all participants. If this is not the case in your study, you need to adjust the code. 

We included data from this study in our pilot data: 

Lawson, R. P., Bisby, J., Nord, C. L., Burgess, N., & Rees, G. (2021). The computational, pharmacological, and physiological determinants of sensory learning under uncertainty. Current Biology, 31(1), 163-172.

Data was downloaded from: https://github.com/BeckyLawson/Propranolol

----------------
RUN THE ANALYSIS
----------------

To run the analysis:
- Make sure that you are inside the directory of the EMBA_PAL repository.
- Specify the model names in the workflow/getModelNames.m
  - Any observation or perception model not contained in the tapas toolbox has to be defined in the 'models' folder
  - All models in the 'models' folder need to have unique names not used in the tapas toolbox
- Specify the necessary parameters in "main.m" and press "run"
  - Do not forget to specify the output directory!
  - Tapas will ask you for an SPM path which is not necessary for this pipeline. Just press ENTER.

----------------
PLOTTING
----------------

All the plotting outside of the main.m workflow is hardcoded and specific to the models that we have used in the EMBA project. 

-------------------
LICENSE INFORMATION
-------------------
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

-------------------
CONTACT INFORMATION
-------------------

If you have any questions or suggestions, please contact the project lead 
Irene Sophia Plank via email at irene.plank@med.uni-muenchen.de

----------
REFERENCES
----------

Please cite the following references when using the EMBA_HGF toolbox:

<INSERT REFERENCE>

Hess, A.J., Iglesias, S., Köchli, L., Marino, S., Müller-Schrader, 
M., Rigoux, L., Mathys, C., Harrison, O.K., Heinzle, J., Frässle, S., 
Stephan, K.E., 2024. Bayesian Workflow for Generative Modeling in Computational Psychiatry. 
*bioRxiv*, https://doi.org/10.1101/2024.02.19.581001

# EMBA: using the Bayesian Brain to investigate the specificity of emotion recognition differences in autism and ADHD

Mutual social interactions require people to be aware of the affective state of their counterpart. An important source for this is facial expressions, which can be indicative of the emotions experienced by the other person. Individuals with autism spectrum disorder (ASD) often have difficulties with this function. Despite the extensive documentation of such differences, it is still unclear which processes underlie attenuated emotion recognition in ASD. In this project, we aim to use a prominent human brain theory called the Bayesian Brain to evaluate the impact of three mechanisms on emotion recognition in individuals with ASD at the neural and behavioural levels: (1) emotional face processing, (2) learning of associations between contextual cues and facial expressions associated with emotions, and (3) biased attention for faces. We also plan to include individuals with attention deficit hyperactivity disorder (ADHD) as clinical controls in addition to a sample of people with no neurodevelopmental disorder. This allows us to determine whether differences in emotion recognition can be attributed to attentional deficits or are unspecific for the included developmental disorders. The results of this project will not only shed more light on the causes of deficits in emotion recognition in people with ASD, but also provide the basis for developing a model of the similarities and differences in processes of the Bayesian Brain in neurodevelopmental disorders.

## Probabilistic associative learning (PAL)

In this repository, we will focus on the probabilistic associative learning (PAL) paradigm. This paradigm was adapted from Lawson et al. (2017, 2021). Contrary to the original task, the participant here is presented with facial expressions associated with fear and happiness and their task is to judge as fast and accurately as possible whether the facial expression is associated with a positively or negatively valenced emotion. Before each face, a high or low tone is played, which can be predictive for the outcome of the following decision. Lawson and colleagues have shown that adults with autism show increased updates of learning rates on the level of the environmental volatility while showing decreased updates of learning rates on the level of the probabilistic outcomes.

Participants also perform three additional paradigms: a dot-probe task to measure face attention bias (FAB), a probabilistic learning paradigm (PAL) and a visual mismatch task (VMM). The preregistrations for this project are on [OSF](https://osf.io/znrht). The preregistrations will be made public when manuscripts or preprints are submitted. 

This repository is a work in progress. The scripts are continuously augmented.

## How to run this analysis

This repository includes scripts for the presentation of the paradigm, preprocessing of the data and analysis. Due to privacy issues, we only share preprocessed and anonymised data. Therefore, we can only guarantee that the final analysis scripts can be run based on this repository: 

* `S1_brms-analyses_PAL.Rmd` : behavioural analysis > run this first
* `S2_brms-analyses_PAL-PUP.Rmd` : pupil size analysis
* `S4_brms-analyses_PAL-HGF.Rmd` : analysis of the Hierarchical Gaussian Filter (HGF) parameters

These scripts may use scripts from the `helpers` folder. There are some absolute paths in these scripts within if statements. Downloading everything in this repository should ensure that these are not executed. 

We also share the models and the results of the simulation-based calibration. **Rerunning these, especially the SBC, can take days depending on the specific model.** Runtime of the scripts using the models and SBC shared in this repository should only take a few minutes. The scripts will create all relevant output that was used in the manuscript. If you need access to other data associated with this project or want to use the stimuli / paradigm, please contact the project lead (Irene Sophia Plank, 10planki@gmail.com). 

### Versions and installation

Each pdf file contains an output of the versions used to run that particular script. It is important to install all packages mentioned in the file before running a specific analysis file. Not all packages can be installed with `install.packages`, please consult the respective installation pages of the packages for more information. If the models are rerun, ensure a valid cmdstanr installation. 

To render the RMarkdown file as a PDF, an installation of pdflatex is mandatory. 

In this project, we used Python 3, MATLAB R2023a and R 4.4.3. 

## Folder content

* `_brms_models` : rds files of the final brms models and the ANOVA to be read into R
* `_brms_SBC_cache` : results of the simulation-based calibration
* `_brms_sens_cache` : results of the sensitivity analysis for the Bayes Factors
* `data` : preprocessed and anonymised data as well as the results from the cross-validation to determine the samples used for the pupil size analysis
* `EMBA_HGF` : Bayesian workflow for the HGF - these scripts are used by `S3_EMBA_HGF_main.m`
* `experiment` : scripts needed to present the experiment as well as the RMarkdown containing all information regarding the stimulus evaluation and selection
* `helper` : scripts to run the simulation-based calibration via bash loops and the script containing the functions for the sensitivity analysis
* `HGF_results` : results of the Bayesian workflow for HGF
* `plots` : some of the plots created for the manuscript
* `prepro` : scripts to preprocess the pupil sizes, the behavioural responses and the demographic as well as questionnaire information collected with CentraXX

## Project members

* Project lead: Irene Sophia Plank
* NEVIA lab PI: Christine M. Falter-Wagner
* Project members (alphabetically): Krasniqi, Kaltrina; Nowak, Julia; Pior, Alexandra; Yurova, Anna; Shi, Zhuanghua

## Licensing

GNU GENERAL PUBLIC LICENSE

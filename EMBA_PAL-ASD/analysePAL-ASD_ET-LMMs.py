# -*- coding: utf-8 -*-
"""
Created on Sun Jan 19 14:00:37 2025

@author: Irene Sophia Plank, 10planki@gmail.com
"""

import glob
import os
from datamatrix import io, DataMatrix, MultiDimensionalColumn, convert
import statsmodels.formula.api as smf
import pandas as pd

####### USE LMS TO DETERMINE THE BETAS FOR EACH PARTICIPANT PER TIMEPOINT

# load the relevant data
path     = '/Users/vilya/.mounty/iplank-ssd/EMBA'
path_out = '/home/emba/Documents/EMBA/PAL-ET_preprocessed/'
path_out = path + '/PAL-ET_preprocessed'
dm = io.readpickle(path_out + '/PAL-ET.pkl')
#df_trj = pd.read_csv('/home/emba/Documents/EMBA/EMBA_PAL_scripts/HGF_results/main/eHGF-C_traj.csv')
df_trj = pd.read_csv(path + '/EMBA_PAL_scripts/EMBA_PAL-ASD/HGF_results/main/HGF-L17_traj.csv')
df_trl = pd.read_csv(path + '/EMBA_PAL_scripts/data/PAL_scheme.csv')

# merge HGF trial output with trial information
df_trl = pd.merge(df_trj, df_trl, how="left")

# get rid of unnecessary columns and rows
del dm.artft, dm.baseline, dm.ds_pupil, dm.key, dm.miss, dm.pupil, dm.z_baseline
dm = dm.exp == {"expected", "unexpected"}
dm = dm.rts > 0
dm = dm.diagnosis == {'COMP', 'ASD'}

# set the formula
lm_form = 'rel_pupil ~ C(exp, Sum) + eps3 + C(emo, Sum) + C(difficulty, Sum) + rts'

# trial duration
dur = dm.bsc_pupil.shape[1]

# initialise the datamatrix
subIDs = list(set(dm.subID))
dm_lms = DataMatrix(length=len(subIDs))
dm_lms.subID   = subIDs
dm_lms.b_Int   = MultiDimensionalColumn(shape=(dur,)) 
dm_lms.b_exp   = MultiDimensionalColumn(shape=(dur,)) 
dm_lms.b_emo   = MultiDimensionalColumn(shape=(dur,)) 
dm_lms.b_diff1 = MultiDimensionalColumn(shape=(dur,)) 
dm_lms.b_diff2 = MultiDimensionalColumn(shape=(dur,))
dm_lms.b_eps3  = MultiDimensionalColumn(shape=(dur,))
dm_lms.b_rts   = MultiDimensionalColumn(shape=(dur,)) 
group = []
i = 0

# loop through the participants
for subID in subIDs:
    
    # track progress
    print('Subject {} of {}'.format(i, dm_lms.shape[0]))

    # extract this subjects data and initialise lists
    dm_subset = dm.subID == subID
    params    = []
    group.append(dm_subset[0].diagnosis)
    
    # loop through the time bins
    for t in range(dur):
        dm_subset.rel_pupil = dm_subset.bsc_pupil[:, t]
        
        # convert to dataframe and merge
        #df_sub = convert.to_pandas(dm_subset)
        d = {}
        for colname, col in dm_subset.columns:
            d[colname] = list(col)
        df_sub = pd.DataFrame(d)
        df_sub = pd.merge(df_sub, df_trl, how="left")
        
        # run the linear model
        lm = smf.ols(formula=lm_form,
                     data=df_sub).fit()
        
        # log all the data
        params.append(lm.params)
    
    # add this subject betas to the datamatrix
    params = pd.DataFrame(params)
    dm_lms[i].b_Int   = params['Intercept']
    dm_lms[i].b_exp   = params['C(exp, Sum)[S.expected]'] 
    dm_lms[i].b_emo   = params['C(emo, Sum)[S.negative]']
    dm_lms[i].b_diff1 = params['C(difficulty, Sum)[S.difficult]']
    dm_lms[i].b_diff2 = params['C(difficulty, Sum)[S.easy]']
    dm_lms[i].b_eps3  = params['eps3']
    dm_lms[i].b_rts   = params['rts']
    
    i = i + 1
    
dm_lms.diagnosis = group

# summarise the timewindows: 1000 - 1250 for expectancy, 500 - 1500 for eps3
# whole thing is 200 samples = 2000ms including 200ms baseline
dm_lms.tiwo_eps3 = dm_lms.b_eps3[:, 70:170]
dm_lms.tiwo_exp  = dm_lms.b_exp[:, 120:145]

# take the mean
dm_lms.tiwo_eps3 = srs.reduce(dm_lms.tiwo_eps3)
dm_lms.tiwo_exp  = srs.reduce(dm_lms.tiwo_exp)

# save all the lm data
io.writepickle(dm_lms, path + '/EMBA_PAL_scripts/data/PAL-ASD-ET_lm-betas.pkl')

# save part of it as csv for further analysis
dm_small = dm_lms[dm_lms.subID, dm_lms.diagnosis, dm_lms.tiwo_eps3, dm_lms.tiwo_exp]
del dm.artft, dm.baseline, dm.ds_pupil, dm.key, dm.miss, dm.pupil, dm.z_baseline

io.writetxt(dm_small, path_out + '/PAL-ASD-ET_lm-betas.csv')

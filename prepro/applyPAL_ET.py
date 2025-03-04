# -*- coding: utf-8 -*-
"""
Created on Sun Jan 19 14:00:37 2025

@author: Irene Sophia Plank, 10planki@gmail.com
"""

import glob
import os
import simple_colors
from datamatrix import io, DataMatrix, MultiDimensionalColumn, convert
import time_series_test_ISP as tst
#import numpy as np
import pandas as pd
from preprocessPAL_ET import PAL_pupilpreprocess
#from datamatrix import operations as ops

####### INITIALISE STUFF

# set the paths
path_ins  = ['/home/emba/Documents/EMBA/BVET/', '/home/emba/Documents/EMBA/BVET-explo/']
path_out  = '/home/emba/Documents/EMBA/PAL-ET_preprocessed/'
path_data = '/home/emba/Documents/EMBA/EMBA_PAL_scripts/data/'

# timing parameters
bsms  = 200  # baseline in ms
trlms = 1800 # trial duraion in ms (face + response window)

# initialise empty datamatrix
dm = DataMatrix(length=0)

####### PREPROCESS PARTICIPANT DAYA ONE BY ONE

for path_in in path_ins:
    # get the filenames
    ls_files = glob.glob(path_in + 'PAL-ET-*.csv')
    # loop through the files
    for et_file in ls_files:
        # extract filename
        file   = os.path.split(et_file)[1].split(sep='PAL-ET-')[1].split(sep='.csv')[0]
        # check if the behavioural exists before preprocessing
        beh_file = glob.glob(path_in + 'PAL-BV-' + file[:10] + '*.csv')
        # check if the subID exists in the datamatrix already [!ISP: MISSING!]
        if len(beh_file) > 0:
            try: 
                dm_sub = PAL_pupilpreprocess(et_file, beh_file[0], 
                                             path_out, path_data,
                                             log=True, lp_filter=4) # Kret: 4 Hz
                # if preprocessing was stopped, None is returned
                if not(dm_sub is None):
                    dm = dm << dm_sub
            except Exception as error:
                print(simple_colors.red("Didn't work': " + file[:10], 'bold'), error)
        else:
            print(simple_colors.red("No behavioural data: " + file[:10], 'bold'))

####### ANONYMISE DATA AND ADD DIAGNOSTIC GROUP

# read the key for anonymisation
df_ano = pd.read_csv(path_ins[0] + '/PID_anonymisation.csv')

# convert the PID column to dataframe and merge > pilots set to NaN
dm_PID = DataMatrix(length=len(dm))
dm_PID.PID = dm.PID
df = convert.to_pandas(dm_PID)
df = pd.merge(df_ano, df, how="right")

# add the information back into the datamatrix
dm.subID     = df['subID']
dm.diagnosis = df['diagnosis']

# get rid of the pseudo-code
del dm.PID

# filter out pilots
dm = dm.subID > 0

# save all preprocessed data in one file
io.writepickle(dm, path_out + 'PAL-ET.pkl')

# read it in again
dm = io.readpickle(path_out + 'PAL-ET.pkl')

####### SUMMARISE DATA PER EXP FOR EACH SUBJECT

# get all the subject IDs
subIDs = set(dm.subID)

# initialise the datamatrix
dm_sum = DataMatrix(length=len(subIDs)*3)
dm_sum.bsc_pupil = MultiDimensionalColumn(shape=(dm.bsc_pupil.shape[1],))

# loop through the subjects
i = 0 # counter index
for sub in subIDs:
    print(sub)
    # loop through conditions
    for level in set(dm.exp):
        # focus on this one subject and condition
        dm_tmp = (dm.subID == sub) & (dm.exp == level)
        # add the information 
        dm_sum[i].subID     = sub
        dm_sum[i].diagnosis = dm_tmp[0].diagnosis
        dm_sum[i].artft     = dm_tmp.artft.mean
        dm_sum[i].miss      = dm_tmp.miss.mean
        dm_sum[i].bsc_pupil = dm_tmp.bsc_pupil.mean
        dm_sum[i].no_trls   = dm_tmp.shape[0]
        dm_sum[i].exp       = level
        i = i + 1

io.writepickle(dm_sum, path_out + 'PAL-ET_sum.pkl')
df = convert.to_pandas(dm_sum)
df.to_csv(path_out + 'PAL-ET_sum.csv')

####### USE CROSS-VALIDATION TO DETERMINE THE DATA FOR THE LMMs

# get rid of the neutral trials and trials without a response
dm_subset = dm.exp == {"expected", "unexpected"}
dm_subset = dm_subset.rts > 0

# get rid of the baseline (200ms = 20 samples) and the pupil reflex to the faces
sr  = 100
ign = int((bsms+800)/(1000/sr))
dm_subset.rel_pupil = dm_subset.bsc_pupil[:, ign:]

# rename to expected
dm_subset.expected = dm_subset.exp

# only keep relevant columns
del dm_subset.artft, dm_subset.baseline, dm_subset.ds_pupil, dm_subset.trl, dm_subset.bsc_pupil
del dm_subset.key, dm_subset.miss, dm_subset.pupil, dm_subset.z_baseline, dm_subset.exp
dm_subset.column_names

results = tst.lmer_crossvalidation_test(dm_subset, split=6,
                                        formula='rel_pupil ~ expected + rts', 
                                        groups='subID', winlen=10, 
                                        con_formula='rel_pupil ~ C(expected, Sum) + rts',
                                        out=path_data + 'CV_pup_sum'
                                        )

print(results)

tst.save(results, path_data + 'CV_pup_sum/results.pkl')

# plot the results with the included function

results = tst.load(path_data + 'CV_pup_sum/results.pkl')

tst.plot(dm_subset, dv='rel_pupil', hue_factor='expected', 
         sampling_freq=100, results=results, x0 = ign/sr)


####### USE LMS TO DETERMINE THE BETAS FOR EACH PARTICIPANT PER TIMEPOINT

# [!ISP]: MISSING



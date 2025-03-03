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
import numpy as np
import pandas as pd
from preprocessPAL_ET import PAL_pupilpreprocess
from matplotlib import pyplot as plt
from datamatrix import operations as ops

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

####### DO SOME PLOTTING

# plot difference between expected and unexpected trials
def plot_series(x, s, color, label):

    se = s.std / np.sqrt(len(s))
    plt.fill_between(x, s.mean-se, s.mean+se, color=color, alpha=.25)
    plt.plot(x, s.mean, color=color, label=label)

x = np.linspace(-bsms, trlms, len(dm[0].bsc_pupil))
dm_ex, dm_un = ops.split(dm_sum.exp, "expected", "unexpected")

plt.figure()
plot_series(x, dm_ex.bsc_pupil, color='green', label='expected (N=%d)' % len(dm_ex))
plot_series(x, dm_un.bsc_pupil, color='blue', label='unexpected (N=%d)' % len(dm_un))
plt.axvline(0, linestyle=':', color='black')
plt.ylabel('Pupil size')
plt.xlabel('Time relative to onset of the face (ms)')
plt.legend(frameon=False, title='')
plt.show()

# average number of artefact samples per condition
data = [dm_ex.artft, dm_un.artft]
fig = plt.figure(figsize =(10, 7))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xticklabels(['expected', 'unexpected'])
bp = ax.boxplot(data)
plt.show()

# average number of missing after correction per condition
data = [dm_ex.miss, dm_un.miss]
fig = plt.figure(figsize =(10, 7))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xticklabels(['expected', 'unexpected'])
bp = ax.boxplot(data)
plt.show()

####### USE CROSS-VALIDATION TO DETERMINE THE DATA FOR THE LMM

# get rid of the neutral trials and trials without a response
dm_subset = dm.exp == {"expected", "unexpected"}
dm_subset = dm_subset.rts > 0

# get rid of the baseline (200ms = 20 samples)
dm_subset.rel_pupil = dm_subset.bsc_pupil[:, 20:]

# only keep relevant columns
del dm_subset.artft, dm_subset.baseline, dm_subset.ds_pupil, dm_subset.bsc_pupil
del dm_subset.key, dm_subset.miss, dm_subset.pupil, dm_subset.z_baseline
dm_subset.column_names

# run the cross validation which saves the relevant csv in the data folder
results = tst.lmer_crossvalidation_test(dm_subset, 
                                        formula='rel_pupil ~ diagnosis * exp + rts', 
                                        groups = 'subID', winlen=1, 
                                        con_formula = 'rel_pupil ~ C(diagnosis, Sum) * C(exp, Sum) + rts',
                                        out = path_data + 'CV_pupilsize'
                                        )

print(results)

####### USE LMS TO DETERMINE THE BETAS FOR EACH PARTICIPANT PER TIMEPOINT

# [!ISP]: MISSING



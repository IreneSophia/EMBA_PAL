# -*- coding: utf-8 -*-
"""
Created on Sun Jan 19 14:00:37 2025

@author: Irene Sophia Plank, 10planki@gmail.com
"""

import glob
import os
import simple_colors
from datamatrix import io, DataMatrix
import time_series_test_ISP as tst
from preprocessPAL_ET import PAL_pupilpreprocess

# set the paths
path_in  = '/home/emba/Documents/EMBA/BVET/'
path_out = '/home/emba/Documents/EMBA/PAL-ET_preprocessed_new/'

# get the filenames
ls_files = glob.glob(path_in + 'PAL-ET-*.csv')
ls_files = ls_files[:10]

# initialise an empty datamatrix
dm = DataMatrix(length=0)

# preprocess participant data one by one
for et_file in ls_files:
    # extract filename
    file   = os.path.split(et_file)[1].split(sep='PAL-ET-')[1].split(sep='.csv')[0]
    # check if the behavioural exists before preprocessing
    beh_file = glob.glob(path_in + 'PAL-BV-' + file[:10] + '*.csv')
    if len(beh_file) > 0:
        try: 
            dm_sub = PAL_pupilpreprocess(et_file, beh_file[0],
                                         out=path_out, lp_filter=None) 
            #dm = dm << dm_sub
        except Exception as error:
            print(simple_colors.red("Didn't work': " + file[:10], 'bold'), error)
    else:
        print(simple_colors.red("No behavioural data: " + file[:10], 'bold'))
        
io.writepickle(dm, path_out + 'PAL-ET.pkl')

# [!ISP]: change rts to 0 if no key was logged


# [!ISP]: summarise for each participant one timeseries per expectancy


# [!ISP]: plot average number of artefact samples per group and condition


####### USE CROSS-VALIDATION TO DETERMINE THE DATA FOR THE LMM

# [!ISP]: only keep relevant columns


# [!ISP]: add an input to the function to determine the output folder
results = tst.lmer_crossvalidation_test(dm, formula='bsc_pupil ~ group * exp', 
                                        groups = 'subject_nr', winlen=1, 
                                        con_formula = 'bsc_pupil ~ C(group, Sum) * C(exp, Sum)'
                                        )

print(results)

# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 19:29:15 2025

@author: Irene Sophia Plank, 10planki@gmail.com
"""

import os 
import simple_colors
import warnings
import pandas as pd
import numpy as np
import traceback 
from datamatrix import DataMatrix, io, MultiDimensionalColumn
from datamatrix import series as srs, operations as ops #functional as fnc, 
from scipy.interpolate import CubicSpline
from matplotlib import pyplot as plt
from time import time
#from matplotlib import lines

# make a function out of this that is then looped over participants
def PAL_pupilpreprocess(et_file, beh_file, out, df_trl, log=False,
                        plotting=True, margin=10, lp_filter=None,
                        sr=500, bcms=200, trlms=1800):
    
    ##### SET UP SOME STUFF
    warnings.filterwarnings('ignore')
    
    # create the output directory if it does not exist
    if not os.path.exists(out):
       os.makedirs(out)
       
    # open log file if set to true
    if log:
        log_path = out + "log.txt"
        log_file = open(log_path, "a")
        
    ###################### DO EVERYTHING AND CATCH POSSIBLE ERRORS
    
    try: 
    
        ##### LOAD DATA
        
        # log starting time
        t = time()
        
        # load the behavioural data
        df_beh = pd.read_csv(beh_file)
        
        # extract the subID
        subID  = df_beh.loc[0,'subID']
        
        # print start time
        print(simple_colors.cyan("START: ", 'bold') + subID)
        if log:
            print(" ", file=log_file)
            print("------------- START: {}".format(subID), file=log_file)
        
        # check if the preprocessed file already exists, if yes load and return
        if os.path.exists(out + 'PAL-ET-' + subID + '.pkl'):
            return io.readpickle(out + 'PAL-ET-' + subID + '.pkl')
        elif os.path.exists(out + 'NOT_PAL-ET-' + subID + '.pkl'):
            return None
        
        # merge both data frames together
        df_beh = df_beh.merge(df_trl)
        
        # change reaction time to 0 if there was no key press
        df_beh.loc[pd.isnull(df_beh.key), 'rt'] = 0
        
        # load the eye tracking data
        try: 
            df_pup = pd.read_csv(et_file, usecols=[' LeftPupilMajorAxis', ' Comment'])            
        except:
            df_pup = pd.read_csv(et_file, usecols=[' PupilMajorAxis', ' Comment'])
            df_pup[' LeftPupilMajorAxis'] = df_pup[' PupilMajorAxis']
        
        # only keep relevant data
        try: 
            start  = df_pup[df_pup[' Comment'].str.contains('pic_1_')].first_valid_index() - 1000
        except:
            start = 0
        try:
            end    = df_pup[df_pup[' Comment'].str.contains('pic_336_')].first_valid_index() + 1000
        except:
            end    = len(df_pup)-1
        if start < 0:
            start = 0
        if end > (len(df_pup)-1):
            end = len(df_pup)-1
        df_pup = df_pup[start:end]
        
        # make sure that the pupil size is a float > strings as NaNs
        df_pup[' LeftPupilMajorAxis'] = pd.to_numeric(df_pup[' LeftPupilMajorAxis'], 
                                                      errors='coerce')
        
        ##### PERFORM ARTEFACT REJECTION
        
        # print how many 0s at the beginning
        pupil = np.copy(df_pup[' LeftPupilMajorAxis'])
        print(subID, "before artefact rejection:", str(round(np.mean(pupil == 0), 3)))
        if log:
            print("{} before artefact rejection: {}".format(subID, round(np.mean(pupil == 0), 3)), file=log_file)
        
        # replace unrealistic values with 0 
        pupil[(pupil > 8) | (pupil < 2)] = 0
        print(subID, "after excluding >8 & <2:", str(round(np.mean(pupil == 0), 3)))
        if log:
            print("{} after excluding >8 & <2: {}".format(subID, round(np.mean(pupil == 0), 3)), file=log_file)
        
        # convert pupil data to velocities
        pupil_smooth = srs._smooth(pupil, winlen=5)
        vt_pupil = np.diff(pupil_smooth, prepend=pupil_smooth[0])
        
        # exclude extreme velocities based on IQR method
        qs = np.quantile(vt_pupil, [0.25, 0.75])
        IQR = qs[1] - qs[0]
        pupil = np.where((vt_pupil < qs[0] - 3 * IQR) | (vt_pupil > qs[1] + 3 * IQR), 
                         0, pupil)
        print(subID, "after excluding extreme velocities:", str(round(np.mean(pupil == 0), 3)))
        if log:
            print("{} after excluding extreme velocities: {}".format(subID, round(np.mean(pupil == 0), 3)), file=log_file)
        
        # find all the 0s
        idx = np.where(pupil == 0)[0]
        
        # split them into chunks
        ls_idx = np.split(idx, np.where(np.diff(idx) > 1)[0] + 1)
        
        # add padding for the blinks
        for i in ls_idx:
            pupil[(i[0]-margin):(i[-1]+margin)] = 0
            
        # exclude extreme pupil sizes based on IQR method
        qs = np.quantile(pupil, [0.25, 0.75])
        IQR = qs[1] - qs[0]
        pupil = np.where((pupil < qs[0] - 1.5 * IQR) | (pupil > qs[1] + 1.5 * IQR), 
                         0, pupil)
        print(subID, "after excluding extreme pupil sizes:", str(round(np.mean(pupil == 0), 3)))
        if log:
            print("{} after excluding extreme pupil sizes: {}".format(subID, round(np.mean(pupil == 0), 3)), file=log_file)
    
        # track how much data is missing
        df_pup['miss']  = pupil == 0
        
        # check if this person is missing more than 1/2 of the data > exclusion
        if np.mean(df_pup['miss']) > 1/2:
            print(simple_colors.red("Too much missing data. " + 
                                    subID + " is excluded from analysis: " + 
                                    str(round(np.mean(pupil == 0), 3)), 
                                    'bold'))
            if log:
                print("EXCLUDED: subject {} is missing too much data: {}".format(subID, round(np.mean(pupil == 0), 3)), file=log_file)
            io.writepickle(None, out + 'NOT_PAL-ET-' + subID + '.pkl')
            return None
        else:
            print(subID, "after padding with margins:", 
                  str(round(np.mean(pupil == 0), 3)))
            if log:
                print("{} after padding with margins: {}".format(subID, round(np.mean(pupil == 0), 3)), file=log_file)
        
        ##### FILTERING AND SMOOTHING
        
        # track signal before correction
        pupil_raw = np.copy(pupil)
        
        # change all zeros to NAs
        pupil[pupil == 0] = np.nan
        
        # check if and what frequency lowpass filter is to be applied
        if lp_filter != None:
            
            # give some feedback
            print("FILTERING THE DATA NOW")
            if log:
                print("{} data is being filtered.".format(subID), file=log_file)
        
            # find all the NAs
            idx = np.where(np.isnan(pupil))[0]
            
            # apply a lowpass filter > first get rid of nans by filling with near values
            df_temp = pd.DataFrame({'pupil' : pupil})
            df_temp = df_temp.fillna(method='ffill').fillna(method='bfill')
            pupil = srs.filter_lowpass(df_temp['pupil'], freq_max = lp_filter, 
                                       sampling_freq = sr)
            
            # set to nan again
            pupil[idx] = np.nan
        
        # smooth the data
        pupil = srs._smooth(pupil)
         
        ##### PERFORM ARTEFACT CORRECTION
        
        # find all the NAs
        idx = np.where(np.isnan(pupil))[0]
        
        # split them into chunks
        ls_idx = np.split(idx, np.where(np.diff(idx) > 1)[0] + 1)
        
        # loop through them and perform artefact correction for all below a max duration
        max_dur = 250
        for i in ls_idx:
            # try to take at least 10 samples before and after
            dur = i[-1] - i[0] + 2
            # check if the duration is acceptable
            if int(dur) <= max_dur:
                # get some points before and after as information for the interpolation
                points = np.array([i[0] - dur, i[0]-1, i[-1]+1, i[-1] + dur])
                         # np.concatenate((np.arange(i[0] - dur, i[0]-1), 
                         #                 np.arange( i[-1]+1, i[-1] + dur)))
                # remove what is outside of the data range
                points[points < 0] = 0
                points[points >= len(pupil)] = len(pupil)-1
                # check for each point if it correspond to NA in the data
                while np.isnan(pupil[points[1]]) & (pupil[points[1]] > pupil[points[0]]):
                    points[1] = points[1] - 1
                while np.isnan(pupil[points[2]]) & (pupil[points[2]] < pupil[points[3]]):
                    points[2] = points[2] + 1
                x = 0
                while np.isnan(pupil[points[0]]) & (x < 10) & (points[3] > 1):
                    points[0] = points[0] - 1
                    x = x + 1
                x = 0
                while np.isnan(pupil[points[3]]) & (x < 10) & (points[3] < (len(pupil)-1)):
                    points[3] = points[3] + 1
                    x = x + 1
                # check if there are still nans
                if np.isnan(np.sum(pupil[points])):
                    # switch to linear interpolation
                    pupil[i] = np.interp(i, points, pupil[points])
                else:
                    # use cubic interpolation and set up the function
                    cs = CubicSpline(points, pupil[points])
                    # apply the function
                    pupil[i] = cs(np.arange(i[0], i[-1]+1, 1))
                    
        # sometimes this creates new extreme values, which we exclude here again
        qs = np.quantile(pupil, [0.25, 0.75])
        IQR = qs[1] - qs[0]
        pupil = np.where((pupil < qs[0] - 1.5 * IQR) | (pupil > qs[1] + 1.5 * IQR), 
                         np.nan, pupil)
                
        # print how many still missing after artefact correction
        print(subID, "after artefact correction:", 
              str(round(np.mean(np.isnan(pupil)), 3)))
        if log:
            print("{} after artefact correction: {}".format(subID, round(np.mean(pupil == 0), 3)), file=log_file)
        
        # add the data back to dataframe
        df_pup['pupil'] = pupil
        
        ##### PLOT THE CONTINUOUS DATA FOR VISUAL INSPECTION
        
        # save the whole time course in plots for visual inspection
        if plotting:
            nr_subpl = 3
            len_plot = 1000
            nr_plots = 6#math.ceil(len(pupil)/(len_plot))
            for i in range(nr_plots): 
                # if first of the subplots, initialise the plot
                if i%nr_subpl == 0:
                    fig, axs = plt.subplots(3, figsize=(20, 15), constrained_layout=True,
                                            sharex=False, sharey=True)
                    fig.supxlabel("Samples", fontsize=14)
                    fig.supylabel("Pupil size (a.u.)", fontsize=14)
                    s = 0
                # create the subplots
                ax = plt.subplot(3, 1, s+1)
                start = (i*len_plot)     
                end   = (i+1)*len_plot-1
                y     = pupil_raw[start:end]
                x     = np.arange(start,end,1)
                if (len_plot-1) > len(y): 
                    x = np.arange(start,start+len(y),1)
                else: 
                    x = np.arange(start,end,1)
                ax.plot(x, pupil_raw[start:end], color="deepskyblue",
                         label="raw", linestyle="solid")
                ax.plot(x, pupil[start:end], color="red",
                         label="reconstructed", linestyle="dashed")
                ax.set_xlim(start, end)
                s     = s + 1
                # if the last subplot has been created, save the figure
                if i%nr_subpl == (nr_subpl-1):
                    plt.ylim(round(np.nanmin(pupil), 1)-0.5, round(np.nanmax(pupil), 1)+0.5)
                    plt.savefig(out + 'PAL-ET-' + subID + '_timecourse_' + str(i+1).zfill(2) + '.jpg')
                    plt.close()
        
        ##### SPLIT INTO TRIALS
        
        # create a data frame with the trigger for the picture
        df_idx = df_pup[df_pup[' Comment'].str.contains('pic')]
        
        # initialise the empty datamatrix
        bcsamples  = int(bcms*(sr/1000))     # size of baseline correction interval in samples
        trlsamples = int(trlms*(sr/1000))    # size of trial duration in samples
        totsamples = bcsamples + trlsamples  # total samples of interest
        dm = DataMatrix(length=336)          # empty datamatrix
        dm.PID   = df_beh['subID']           # subject ID
        dm.trl   = df_beh['trl']             # trial number
        dm.key   = df_beh['key']             # response key pressed
        dm.rts   = df_beh['rt']              # reaction time
        dm.exp   = df_beh['expected']        # condition of expectancy
        dm.pupil = MultiDimensionalColumn(shape=(totsamples,)) # pupil data 
        dm.artft = 0                         # artefacts
        dm.miss  = 0                         # missing samples after artefact correction
        
        # loop through the trials and build the dataframe by adding information
        for index, row in df_idx.iterrows():
            
            # extract the trial number
            trl = int(row[' Comment'].split('_')[1])
            
            # add it to the datamatrix
            dm.pupil[trl - 1] = df_pup.loc[(index-bcsamples):(index+trlsamples-1),'pupil']
            dm.artft[trl - 1] = np.mean(df_pup.loc[(index-bcsamples):(index+trlsamples-1),'miss'])
            dm.miss[trl - 1]  = np.mean(np.isnan(dm.pupil[trl - 1]))
        
        # perform downsampling
        dm.ds_pupil = srs.downsample(dm.pupil, by=5)
        
        # exclude trials with more than 1/3 of datapoints missing after interpolation
        dm_subset = dm.miss < 1/3
        
        # perform baseline correction
        bcsamples = int(bcsamples/5) # same as downsampling
        dm_subset.bsc_pupil = srs.baseline(
            dm_subset.ds_pupil, dm_subset.ds_pupil, bl_start=0, bl_end=bcsamples)
        
        # get baseline pupil size
        dm_subset.baseline = srs.reduce(dm_subset.ds_pupil[:, 0:bcsamples])
        
        # plot the z-scores of the baseline pupil size
        dm_subset.z_baseline = ops.z(dm_subset.baseline)
        if plotting:
            plt.hist(dm_subset.z_baseline)
            plt.plot((-2, -2), (0, 50), scaley = False)
            plt.plot((+2, +2), (0, 50), scaley = False)
            plt.ylim(0,50)
            plt.savefig(out + 'PAL-ET-' + subID + '_baselinepupils.jpg')
            plt.close()
        
        # exclude outliers based on baseline pupil size
        dm_subset = dm_subset.z_baseline >= -2
        dm_subset = dm_subset.z_baseline <= +2
        
        # plot the traces of baseline corrected trials > any outliers?
        if plotting: 
            fig, axs = plt.subplots(1, 2, figsize=(
                8, 5), constrained_layout=True, sharex=True, sharey=True)
            fig.supxlabel('Time', color='black', fontsize=14)
            fig.supylabel('Pupil size (mm)', color='black', fontsize=14)
            plt.subplot(121)  # original
            plt.title(r"$\bf{" + 'a) ' + "}$" + ' Original', fontsize=14, loc='left')
            plt.plot(dm_subset.pupil.plottable)
            plt.subplot(122)  # baseline-corrected
            plt.title(r"$\bf{" + 'b) ' + "}$" +
                      ' Baseline-corrected', fontsize=14, loc='left')
            plt.plot(dm_subset.bsc_pupil.plottable)
            plt.axvline(bcsamples, linestyle=':', color='black')
            plt.xlim(0, 210)
            plt.ylim(-3, 8)
            plt.savefig(out + 'PAL-ET-' + subID + '_traceplot.jpg')
            plt.close()
        
        # save the participant data 
        if len(dm_subset)/336 < 2/3: 
            print(simple_colors.red("Too many missing trials. " + subID + " is excluded from analysis.", 'bold'))
            io.writepickle(dm_subset, out + 'NOT_PAL-ET-' + subID + '.pkl')
            if log:
                print("EXCLUDED: subject {} is missing too many trials.".format(subID), file=log_file)
            return None
        else:
            io.writepickle(dm_subset, out + 'PAL-ET-' + subID + '.pkl')
            
        # print end time
        elapsed = time() - t
        print(simple_colors.cyan("FINISHED ", 'bold') + "in " + str(round(elapsed)) + " sec.")
        if log:
            print("FINISHED in {} sec.".format(round(elapsed)), file=log_file)
            log_file.close()
            
        return dm_subset
    
    except:
        traceback.print_exc(file=log_file)
        log_file.close()
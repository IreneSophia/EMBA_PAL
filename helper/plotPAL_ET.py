#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  4 08:54:27 2025

Misc plotting for pupil data in the EMBA_PAL projects

@author: Irene Sophia Plank, 10planki@gmail.com
"""

import glob
import os
from datamatrix import io, DataMatrix, MultiDimensionalColumn, convert
import numpy as np
import pandas as pd
#from preprocessPAL_ET import PAL_pupilpreprocess
from matplotlib import pyplot as plt
from datamatrix import operations as ops

####### INITIALISE STUFF

# timing parameters
bsms  = 200  # baseline in ms
trlms = 1800 # trial duraion in ms (face + response window)

# load the datamatrices
dm     = io.readpickle(path_out + '/PAL-ET.pkl')
dm_sum = io.readpickle(path_out + '/PAL-ET_sum.pkl')

####### DO SOME PLOTTING

# function to plot with standard errors
def plot_series(x, s, color, label):
    se = s.std / np.sqrt(len(s))
    plt.fill_between(x, s.mean-se, s.mean+se, color=color, alpha=.25)
    plt.plot(x, s.mean, color=color, label=label)

# time line
x = np.linspace(-bsms, trlms, len(dm[0].bsc_pupil))

# plot absolute values for everything
y  = dm.bsc_pupil.mean
y  = abs(y)
plt.figure()
plt.plot(x, y, color='blue')
plt.axvline(0, linestyle=':', color='black')
plt.ylabel('abs(pupil size)')
plt.xlabel('Time relative to onset of the face (ms)')
plt.legend(frameon=False, title='')
plt.show()

# separate into expected and unexpected trials
dm_ex, dm_un = ops.split(dm_sum.exp, "expected", "unexpected")
dm_ex_ADHD, dm_ex_ASD, dm_ex_BOTH, dm_ex_COMP = ops.split(dm_ex.diagnosis, "ADHD", "ASD", "BOTH", "COMP")
dm_un_ADHD, dm_un_ASD, dm_un_BOTH, dm_un_COMP = ops.split(dm_un.diagnosis, "ADHD", "ASD", "BOTH", "COMP")

# expectancy effect
plt.figure()
plot_series(x, dm_ex.bsc_pupil, color='green', label='expected (N=%d)' % len(dm_ex))
plot_series(x, dm_un.bsc_pupil, color='blue', label='unexpected (N=%d)' % len(dm_un))
plt.axvline(0,    color='black')
plt.ylabel('Pupil size')
plt.xlabel('Time relative to onset of the face (ms)')
plt.legend(frameon=True, title='')
plt.show()

# separate plots per group
plt.figure()
plot_series(x, dm_ex_ADHD.bsc_pupil, color='green', label='ADHD expected (N=%d)' % len(dm_ex_ADHD))
plot_series(x, dm_un_ADHD.bsc_pupil, color='blue', label='ADHD unexpected (N=%d)' % len(dm_un_ADHD))
plt.axvline(0, color='black')
plt.ylabel('Pupil size')
plt.xlabel('Time relative to onset of the face (ms)')
plt.legend(frameon=True, title='')
plt.show()
plt.figure()
plot_series(x, dm_ex_ASD.bsc_pupil, color='green', label='ASD expected (N=%d)' % len(dm_ex_ASD))
plot_series(x, dm_un_ASD.bsc_pupil, color='blue', label='ASD unexpected (N=%d)' % len(dm_un_ASD))
plt.axvline(0, color='black')
plt.ylabel('Pupil size')
plt.xlabel('Time relative to onset of the face (ms)')
plt.legend(frameon=True, title='')
plt.show()
plt.figure()
plot_series(x, dm_ex_BOTH.bsc_pupil, color='green', label='BOTH expected (N=%d)' % len(dm_ex_BOTH))
plot_series(x, dm_un_BOTH.bsc_pupil, color='blue', label='BOTH unexpected (N=%d)' % len(dm_un_BOTH))
plt.axvline(0, color='black')
plt.ylabel('Pupil size')
plt.xlabel('Time relative to onset of the face (ms)')
plt.legend(frameon=True, title='')
plt.show()
plt.figure()
plot_series(x, dm_ex_COMP.bsc_pupil, color='green', label='COMP expected (N=%d)' % len(dm_ex_COMP))
plot_series(x, dm_un_COMP.bsc_pupil, color='blue', label='COMP unexpected (N=%d)' % len(dm_un_COMP))
plt.axvline(0, color='black')
plt.ylabel('Pupil size')
plt.xlabel('Time relative to onset of the face (ms)')
plt.legend(frameon=True, title='')
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

# LM beta results

path_data = path + '/EMBA_PAL_scripts/data'
dm_lms = io.readpickle(path_data + '/PAL-ASD-ET_lm-betas.pkl')
dm_lms = io.readpickle(path_data + '/PAL-ADHD-ET_lm-betas.pkl')

dm_ADHD, dm_ASD, dm_BOTH, dm_COMP = ops.split(dm_lms.diagnosis, "ADHD", "ASD", "BOTH", "COMP")

plt.figure()
#plot_series(x, dm_ASD.b_eps2,  color='green', label='ASD')
plot_series(x, dm_COMP.b_eps2, color='blue', label='COMP')
plot_series(x, dm_ADHD.b_eps2,  color='red', label='ADHD')
plot_series(x, dm_BOTH.b_eps2, color='yellow', label='BOTH')
plt.axvline(0,    color='black')
plt.ylabel('b_eps2')
plt.xlabel('Time relative to onset of the face (ms)')
plt.legend(frameon=False, title='')
plt.show()

plt.figure()
plot_series(x, dm_ASD.b_eps3,  color='green', label='ASD')
plot_series(x, dm_COMP.b_eps3, color='blue', label='COMP')
#plot_series(x, dm_ADHD.b_eps3,  color='red', label='ADHD')
#plot_series(x, dm_BOTH.b_eps3, color='yellow', label='BOTH')
plt.axvline(0,    color='black')
plt.axvline(500, color='black', linestyle=':')
plt.axvline(1500, color='black', linestyle=':')
plt.ylabel('b_eps3')
plt.xlabel('Time relative to onset of the face (ms)')
plt.legend(frameon=False, title='')
plt.show()

plt.figure()
#plot_series(x, dm_ASD.b_exp,  color='green', label='ASD')
plot_series(x, dm_COMP.b_exp, color='blue', label='COMP')
plot_series(x, dm_ADHD.b_exp,  color='red', label='ADHD')
plot_series(x, dm_BOTH.b_exp, color='yellow', label='BOTH')
plt.axvline(0,    color='black')
plt.axvline(1000, color='black', linestyle=':')
plt.axvline(1250, color='black', linestyle=':')
plt.ylabel('b_exp')
plt.xlabel('Time relative to onset of the face (ms)')
plt.legend(frameon=False, title='')
plt.show()

plt.figure()
plot_series(x, dm_COMP.b_mu3, color='blue', label='COMP')
plot_series(x, dm_ADHD.b_mu3,  color='red', label='ADHD')
plot_series(x, dm_BOTH.b_mu3, color='yellow', label='BOTH')
plt.axvline(0,    color='black')
plt.ylabel('b_mu3')
plt.xlabel('Time relative to onset of the face (ms)')
plt.legend(frameon=False, title='')
plt.show()

plt.figure()
plot_series(x, dm_lms.b_exp,  color='green', label='all')
plt.axvline(0,    color='black')
plt.axvline(1000, color='black', linestyle=':')
plt.axvline(1250, color='black', linestyle=':')
plt.ylabel('b_exp')
plt.xlabel('Time relative to onset of the face (ms)')
plt.legend(frameon=False, title='')
plt.show()

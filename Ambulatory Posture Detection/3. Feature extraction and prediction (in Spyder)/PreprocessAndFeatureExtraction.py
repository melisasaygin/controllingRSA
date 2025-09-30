# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 13:53:30 2024

@author: sve254
"""
# clear console and environment
#%reset -f 
import os
#os.system('cls')
import gc
#gc.collect()

## install packages needed
# !pip install pandas
# !pip install numpy
# !pip install scipy
# !pip install mne
# !pip install librosa
# !pip install pickle

# load libraries
import pandas as pd # handling dataframes
import numpy as np # array handling
from scipy import signal
from scipy.stats import skew, kurtosis # skewness and kurtosis
import datetime
import statistics
from math import sqrt, acos
import librosa
import h5py

# set directory
os.chdir("C:/Users/msa583/Documents/New Posture Detections - Part 2/82556/day2") #change here
folderpathData = "C:/Users/msa583/Documents/New Posture Detections - Part 2/82556/day2"
folderpathTimestamps = "C:/Users/msa583/Documents/New Posture Detections - Part 2/82556/day2"

def StartEndActivities(activityArray):
    activityArray=np.copy(activityArray).astype(object)
    # find indices where activity changes
    activityStart = np.where(activityArray[:-1] != activityArray[1:])[0] + 1
    # add start and end indices for first and last activity
    activityStart = np.concatenate(([0], activityStart, [len(activityArray)]))
    # get run lengths of each activity
    activityLength = np.diff(activityStart)
    # activity names
    activityName = activityArray[activityStart[:-1]]
    # create output array
    activitiesOutput = np.zeros((len(activityLength), 4), dtype=object)
    # populate output array with start and end indices and activity names
    for i in range(len(activityLength)):
        activitiesOutput[i, 0] = activityStart[i]
        activitiesOutput[i, 1] = activityStart[i] + activityLength[i] - 1
        activitiesOutput[i, 2] = activityLength[i]
        activitiesOutput[i, 3] = activityName[i]    
    return activitiesOutput

# preprocess raw accelerometer data and feature extraction ----

# set path, find files and extract participant numbers
fsOld = 1000
fs = 100
filepaths = [os.path.join(folderpathData, f) for f in os.listdir(folderpathData) if f.endswith('.h5')]
files = [f for f in os.listdir(folderpathData) if f.endswith('.h5')]
participants = [f[0:5] for f in files] # maybe adjust, depending on the naming of your csv-files

# set epoch length , set window length to 60 sec (for minute-by-minute prediction)
epochLength = 6
nSamples = int(epochLength*fs)
overlap = 0
windowLength = 60
epochsPerWindow = int((windowLength/epochLength) * (1/(1-overlap)))

# filter characteristics (2nd order 0.1 Hz high-pass butterworth)
orderGravity = 1
cutoffFreqGravity = 0.1
criticalFreqGravity = cutoffFreqGravity/(fsOld/2) # cutoff frequency divided by the nyquist frequency
bGravity, aGravity = signal.butter(orderGravity, criticalFreqGravity, 'high')

# filter characteristics (4th order low-pass butterworth)
order = 2
cutoffFreq = 5
criticalFreq = cutoffFreq/(fsOld/2) # cutoff frequency divided by the nyquist frequency
b, a = signal.butter(order, criticalFreq, 'low')

# read in timestamps file (for later)
# startingTimes = pd.read_csv('StartingTimes.csv', dtype={'Participant': str})
# startingTimes['starting_time_a'] = pd.to_datetime(startingTimes['starting_time_a'], format='%H:%M:%S')
    
# loop over all participants to extract features and create one matrix
for p in range(len(participants)):
    print(participants[p])
    print(p+1, '/', len(participants))
    
    # load accelerometer data  
    with h5py.File(filepaths[p], 'r') as f:
        MXR = pd.Series(f['trimmed_MXR'][:].flatten())
        MYR = pd.Series(f['trimmed_MYR'][:].flatten())
        MZR = pd.Series(f['trimmed_MZR'][:].flatten())
        GyroX = pd.Series(f['trimmed_GyroX'][:].flatten())
        GyroY = pd.Series(f['trimmed_GyroY'][:].flatten())
        GyroZ = pd.Series(f['trimmed_GyroZ'][:].flatten())
    
    # filter acceleration signals with a 2nd order high-pass butterworth
    MXRGravityCorr = signal.filtfilt(bGravity, aGravity, MXR)
    MYRGravityCorr = signal.filtfilt(bGravity, aGravity, MYR)
    MZRGravityCorr = signal.filtfilt(bGravity, aGravity, MZR)
    GyroXGravityCorr = signal.filtfilt(bGravity, aGravity, GyroX)
    GyroYGravityCorr = signal.filtfilt(bGravity, aGravity, GyroY)
    GyroZGravityCorr = signal.filtfilt(bGravity, aGravity, GyroZ)
    
    # filter acceleration signals with a 4th order low-pass butterworth
    MXR = signal.filtfilt(b, a, MXR)
    MYR = signal.filtfilt(b, a, MYR)
    MZR = signal.filtfilt(b, a, MZR)
    GyroX = signal.filtfilt(b, a, GyroX)
    GyroY = signal.filtfilt(b, a, GyroY)
    GyroZ = signal.filtfilt(b, a, GyroZ)
    MXRGravityCorr = signal.filtfilt(b, a, MXRGravityCorr)
    MYRGravityCorr = signal.filtfilt(b, a, MYRGravityCorr)
    MZRGravityCorr = signal.filtfilt(b, a, MZRGravityCorr)
    GyroXGravityCorr = signal.filtfilt(b, a, GyroXGravityCorr)
    GyroYGravityCorr = signal.filtfilt(b, a, GyroYGravityCorr)
    GyroZGravityCorr = signal.filtfilt(b, a, GyroZGravityCorr)
    
    # resample to a lower sampling frequency using 8th order Chebyshev infinite
    # impulse response (IIR) filter
    MXR = signal.decimate(MXR, round(fsOld/fs), ftype='iir', zero_phase=True)
    MYR = signal.decimate(MYR, round(fsOld/fs), ftype='iir', zero_phase=True)
    MZR = signal.decimate(MZR, round(fsOld/fs), ftype='iir', zero_phase=True)
    GyroX = signal.decimate(GyroX, round(fsOld/fs), ftype='iir', zero_phase=True)
    GyroY = signal.decimate(GyroY, round(fsOld/fs), ftype='iir', zero_phase=True)
    GyroZ = signal.decimate(GyroZ, round(fsOld/fs), ftype='iir', zero_phase=True)
    MXRGravityCorr = signal.decimate(MXRGravityCorr, round(fsOld/fs), ftype='iir', zero_phase=True)
    MYRGravityCorr = signal.decimate(MYRGravityCorr, round(fsOld/fs), ftype='iir', zero_phase=True)
    MZRGravityCorr = signal.decimate(MZRGravityCorr, round(fsOld/fs), ftype='iir', zero_phase=True)
    GyroXGravityCorr = signal.decimate(GyroXGravityCorr, round(fsOld/fs), ftype='iir', zero_phase=True)
    GyroYGravityCorr = signal.decimate(GyroYGravityCorr, round(fsOld/fs), ftype='iir', zero_phase=True)
    GyroZGravityCorr = signal.decimate(GyroZGravityCorr, round(fsOld/fs), ftype='iir', zero_phase=True)
    
    # create vector magnitude
    vectorMagnitude = np.sqrt(MXR**2 + MYR**2 + MZR**2)
    vectorMagnitudeGravityCorr = np.sqrt(MXRGravityCorr**2 + MYRGravityCorr**2 + MZRGravityCorr**2)
    vectorMagnitudeGyro = np.sqrt(GyroX**2 + GyroY**2 + GyroZ**2)
    
    # create time vector
    #start_time = startingTimes.loc[startingTimes['Participant'] == participants[p], 'starting_time_a'].iloc[0]
    start_time =  pd.to_datetime('2024-11-16 04:42:48.500', format="%Y-%m-%d %H:%M:%S.%f")
    end_time = start_time + datetime.timedelta(seconds=(len(MXR)-1) * (1/fs))
    time = [start_time + datetime.timedelta(seconds=i*(1/fs)) for i in range(len(MXR))]
    time = [t.strftime("%Y-%m-%d %H:%M:%S.%f") for t in time]
    time = pd.to_datetime(pd.Series(time))
    
    # find periods where data is missing
    indexMissing = np.where((MXR < -3) & (MYR < -3) & (MZR < -3))[0]
    if len(indexMissing)>0:
        start_index = indexMissing[0]
        end_index = indexMissing[0]
        periods = []
        for im in range(1, len(indexMissing)):
            if indexMissing[im] - indexMissing[im-1] > 1:
                # found the end of a period
                end_index = indexMissing[im-1]
                periods.append((start_index, end_index))
                start_index = indexMissing[im]
            else:
                end_index = indexMissing[im]
        # last period
        periods.append((start_index, end_index)) 
        # print number and length of periods
        print('found',len(periods),'missing period(s):')
        for mp in range(len(periods)):
            print((periods[mp][1] - periods[mp][0]) / 100, 'seconds')
        del start_index, end_index, indexMissing
    
    # divide in 6s epochs with 50% overlap
    nEpochs = int(round(len(MXR)/nSamples))
    nEpochsOverlap = int(nEpochs*(1/(1-overlap))-2)
    
    # create empty matrices to fill with features  
    # time
    timeEpochs = np.empty((nEpochsOverlap,1),dtype='datetime64[ns]')
    timeEpochs.fill(np.datetime64('NaT'))
    # features accelerometer
    averagePerAxis = np.empty((nEpochsOverlap,4))
    stdPerAxis = np.empty((nEpochsOverlap,4))
    variancePerAxis = np.empty((nEpochsOverlap,4))
    maxPerAxis = np.empty((nEpochsOverlap,4))
    maxNormPerAxis = np.empty((nEpochsOverlap,4))
    minPerAxis = np.empty((nEpochsOverlap,4))
    minNormPerAxis = np.empty((nEpochsOverlap,4))
    rangePerAxis = np.empty((nEpochsOverlap,4))
    kurtosisPerAxis = np.empty((nEpochsOverlap,4))
    skewnessPerAxis = np.empty((nEpochsOverlap,4))
    corrPerAxis = np.empty((nEpochsOverlap,3))
    inclPerAxis = np.empty((nEpochsOverlap,3))
    zcrPerAxis = np.empty((nEpochsOverlap,3))
    # features gyroscope
    averagePerAxisGyro = np.empty((nEpochsOverlap,4))
    stdPerAxisGyro = np.empty((nEpochsOverlap,4))
    variancePerAxisGyro = np.empty((nEpochsOverlap,4))
    maxPerAxisGyro = np.empty((nEpochsOverlap,4))
    maxNormPerAxisGyro = np.empty((nEpochsOverlap,4))
    minPerAxisGyro = np.empty((nEpochsOverlap,4))
    minNormPerAxisGyro = np.empty((nEpochsOverlap,4))
    rangePerAxisGyro = np.empty((nEpochsOverlap,4))
    kurtosisPerAxisGyro = np.empty((nEpochsOverlap,4))
    skewnessPerAxisGyro = np.empty((nEpochsOverlap,4))
    corrPerAxisGyro = np.empty((nEpochsOverlap,3))
    inclPerAxisGyro = np.empty((nEpochsOverlap,3))
    zcrPerAxisGyro = np.empty((nEpochsOverlap,3))
    # empty list to save all features and processed data
    filteredData = [None] * nEpochsOverlap
    features = [None] * nEpochsOverlap
    #dominantFreqPerAxis = matrix(NA,nEpochsOverlap,3)
    #powerDominantFreqPerAxis = matrix(NA,nEpochsOverlap,3)
    
    # loop over all epochs to extract features and filtered data per epoch
    for j in range(nEpochsOverlap):
        # cut out 6 second epochs of data with 50 % overlap
        startEpoch = int(j * (nSamples - (overlap*nSamples)))
        endEpoch = int(startEpoch + nSamples)
        MXRtemp = MXR[startEpoch:endEpoch]
        MYRtemp = MYR[startEpoch:endEpoch]
        MZRtemp = MZR[startEpoch:endEpoch]
        VMtemp = vectorMagnitudeGravityCorr[startEpoch:endEpoch]
        GyroXtemp = GyroX[startEpoch:endEpoch]
        GyroYtemp = GyroY[startEpoch:endEpoch]
        GyroZtemp = GyroZ[startEpoch:endEpoch]
        GyroVMtemp = vectorMagnitudeGyro[startEpoch:endEpoch]
        # time
        timeEpochs[j] = time[startEpoch]
        # mean
        averagePerAxis[j,0] = statistics.mean(MXRtemp)
        averagePerAxis[j,1] = statistics.mean(MYRtemp)
        averagePerAxis[j,2] = statistics.mean(MZRtemp)
        averagePerAxis[j,3] = statistics.mean(VMtemp)
        # standard deviation
        stdPerAxis[j,0] = statistics.stdev(MXRtemp)
        stdPerAxis[j,1] = statistics.stdev(MYRtemp)
        stdPerAxis[j,2] = statistics.stdev(MZRtemp)
        stdPerAxis[j,3] = statistics.stdev(VMtemp)
        # variance
        variancePerAxis[j,0] = statistics.variance(MXRtemp)
        variancePerAxis[j,1] = statistics.variance(MYRtemp)
        variancePerAxis[j,2] = statistics.variance(MZRtemp)
        variancePerAxis[j,3] = statistics.variance(VMtemp)
        # maximum
        maxPerAxis[j,0] = max(MXRtemp)
        maxPerAxis[j,1] = max(MYRtemp)
        maxPerAxis[j,2] = max(MZRtemp)
        maxPerAxis[j,3] = max(VMtemp)
        # maximum normalized
        maxNormPerAxis[j,0] = max(MXRtemp - statistics.mean(MXRtemp))
        maxNormPerAxis[j,1] = max(MYRtemp - statistics.mean(MYRtemp))
        maxNormPerAxis[j,2] = max(MZRtemp - statistics.mean(MZRtemp))
        maxNormPerAxis[j,3] = max(VMtemp - statistics.mean(VMtemp))
        # minimum
        minPerAxis[j,0] = min(MXRtemp)
        minPerAxis[j,1] = min(MYRtemp)
        minPerAxis[j,2] = min(MZRtemp)
        minPerAxis[j,3] = min(VMtemp)
        # maximum normalized
        minNormPerAxis[j,0] = min(MXRtemp - statistics.mean(MXRtemp))
        minNormPerAxis[j,1] = min(MYRtemp - statistics.mean(MYRtemp))
        minNormPerAxis[j,2] = min(MZRtemp - statistics.mean(MZRtemp))
        minNormPerAxis[j,3] = min(VMtemp - statistics.mean(VMtemp))
        # range
        rangePerAxis[j,0] = maxPerAxis[j,0] - minPerAxis[j,0]
        rangePerAxis[j,1] = maxPerAxis[j,1] - minPerAxis[j,1]
        rangePerAxis[j,2] = maxPerAxis[j,2] - minPerAxis[j,2]
        rangePerAxis[j,3] = maxPerAxis[j,3] - minPerAxis[j,3]
        # kurtosis
        kurtosisPerAxis[j,0] = kurtosis(MXRtemp)
        kurtosisPerAxis[j,1] = kurtosis(MYRtemp)
        kurtosisPerAxis[j,2] = kurtosis(MZRtemp)
        kurtosisPerAxis[j,3] = kurtosis(VMtemp)
        # skewness
        skewnessPerAxis[j,0] = skew(MXRtemp)
        skewnessPerAxis[j,1] = skew(MYRtemp)
        skewnessPerAxis[j,2] = skew(MZRtemp)
        skewnessPerAxis[j,3] = skew(VMtemp)
        # pairwise correlation
        corrPerAxis[j,0] = np.corrcoef(MXRtemp,MYRtemp)[0,1] # pairwise X-Y
        corrPerAxis[j,1] = np.corrcoef(MXRtemp,MZRtemp)[0,1] # pairwise X-Z
        corrPerAxis[j,2] = np.corrcoef(MYRtemp,MZRtemp)[0,1] # pairwise Y-Z
        # inclination per axis
        inclPerAxis[j,0] = np.rad2deg(acos(averagePerAxis[j,0]/sqrt(averagePerAxis[j,0]**2+averagePerAxis[j,1]**2+averagePerAxis[j,2]**2)))
        inclPerAxis[j,1] = np.rad2deg(acos(averagePerAxis[j,1]/sqrt(averagePerAxis[j,0]**2+averagePerAxis[j,1]**2+averagePerAxis[j,2]**2)))
        inclPerAxis[j,2] = np.rad2deg(acos(averagePerAxis[j,2]/sqrt(averagePerAxis[j,0]**2+averagePerAxis[j,1]**2+averagePerAxis[j,2]**2)))
        ## zero crossing rate, normalize so average value is zero
        zcrPerAxis[j,0] = librosa.feature.zero_crossing_rate(MXRtemp - statistics.mean(MXRtemp))[0,1]
        zcrPerAxis[j,1] = librosa.feature.zero_crossing_rate(MYRtemp - statistics.mean(MYRtemp))[0,1]
        zcrPerAxis[j,2] = librosa.feature.zero_crossing_rate(MZRtemp - statistics.mean(MZRtemp))[0,1]
        
        # mean
        averagePerAxisGyro[j,0] = statistics.mean(GyroXtemp)
        averagePerAxisGyro[j,1] = statistics.mean(GyroYtemp)
        averagePerAxisGyro[j,2] = statistics.mean(GyroZtemp)
        averagePerAxisGyro[j,3] = statistics.mean(GyroVMtemp)
        # standard deviation
        stdPerAxisGyro[j,0] = statistics.stdev(GyroXtemp)
        stdPerAxisGyro[j,1] = statistics.stdev(GyroYtemp)
        stdPerAxisGyro[j,2] = statistics.stdev(GyroZtemp)
        stdPerAxisGyro[j,3] = statistics.stdev(GyroVMtemp)
        # variance
        variancePerAxisGyro[j,0] = statistics.variance(GyroXtemp)
        variancePerAxisGyro[j,1] = statistics.variance(GyroYtemp)
        variancePerAxisGyro[j,2] = statistics.variance(GyroZtemp)
        variancePerAxisGyro[j,3] = statistics.variance(GyroVMtemp)
        # maximum
        maxPerAxisGyro[j,0] = max(GyroXtemp)
        maxPerAxisGyro[j,1] = max(GyroYtemp)
        maxPerAxisGyro[j,2] = max(GyroZtemp)
        maxPerAxisGyro[j,3] = max(GyroVMtemp)
        # maximum normalized
        maxNormPerAxisGyro[j,0] = max(GyroXtemp - statistics.mean(GyroXtemp))
        maxNormPerAxisGyro[j,1] = max(GyroYtemp - statistics.mean(GyroYtemp))
        maxNormPerAxisGyro[j,2] = max(GyroZtemp - statistics.mean(GyroZtemp))
        maxNormPerAxisGyro[j,3] = max(GyroVMtemp - statistics.mean(GyroVMtemp))
        # minimum
        minPerAxisGyro[j,0] = min(GyroXtemp)
        minPerAxisGyro[j,1] = min(GyroYtemp)
        minPerAxisGyro[j,2] = min(GyroZtemp)
        minPerAxisGyro[j,3] = min(GyroVMtemp)
        # maximum normalized
        minNormPerAxisGyro[j,0] = min(GyroXtemp - statistics.mean(GyroXtemp))
        minNormPerAxisGyro[j,1] = min(GyroYtemp - statistics.mean(GyroYtemp))
        minNormPerAxisGyro[j,2] = min(GyroZtemp - statistics.mean(GyroZtemp))
        minNormPerAxisGyro[j,3] = min(GyroVMtemp - statistics.mean(GyroVMtemp))
        # range
        rangePerAxisGyro[j,0] = maxPerAxisGyro[j,0] - minPerAxisGyro[j,0]
        rangePerAxisGyro[j,1] = maxPerAxisGyro[j,1] - minPerAxisGyro[j,1]
        rangePerAxisGyro[j,2] = maxPerAxisGyro[j,2] - minPerAxisGyro[j,2]
        rangePerAxisGyro[j,3] = maxPerAxisGyro[j,3] - minPerAxisGyro[j,3]
        # kurtosis
        kurtosisPerAxisGyro[j,0] = kurtosis(GyroXtemp)
        kurtosisPerAxisGyro[j,1] = kurtosis(GyroYtemp)
        kurtosisPerAxisGyro[j,2] = kurtosis(GyroZtemp)
        kurtosisPerAxisGyro[j,3] = kurtosis(GyroVMtemp)
        # skewness
        skewnessPerAxisGyro[j,0] = skew(GyroXtemp)
        skewnessPerAxisGyro[j,1] = skew(GyroYtemp)
        skewnessPerAxisGyro[j,2] = skew(GyroZtemp)
        skewnessPerAxisGyro[j,3] = skew(GyroVMtemp)
        # pairwise correlation
        corrPerAxisGyro[j,0] = np.corrcoef(GyroXtemp,GyroYtemp)[0,1] # pairwise X-Y
        corrPerAxisGyro[j,1] = np.corrcoef(GyroXtemp,GyroZtemp)[0,1] # pairwise X-Z
        corrPerAxisGyro[j,2] = np.corrcoef(GyroYtemp,GyroZtemp)[0,1] # pairwise Y-Z
        # inclination per axis
        inclPerAxisGyro[j,0] = np.rad2deg(acos(averagePerAxisGyro[j,0]/sqrt(averagePerAxisGyro[j,0]**2+averagePerAxisGyro[j,1]**2+averagePerAxisGyro[j,2]**2)))
        inclPerAxisGyro[j,1] = np.rad2deg(acos(averagePerAxisGyro[j,1]/sqrt(averagePerAxisGyro[j,0]**2+averagePerAxisGyro[j,1]**2+averagePerAxisGyro[j,2]**2)))
        inclPerAxisGyro[j,2] = np.rad2deg(acos(averagePerAxisGyro[j,2]/sqrt(averagePerAxisGyro[j,0]**2+averagePerAxisGyro[j,1]**2+averagePerAxisGyro[j,2]**2)))
        ## zero crossing rate, normalize so average value is zero
        zcrPerAxisGyro[j,0] = librosa.feature.zero_crossing_rate(GyroXtemp - statistics.mean(GyroXtemp))[0,1]
        zcrPerAxisGyro[j,1] = librosa.feature.zero_crossing_rate(GyroYtemp - statistics.mean(GyroYtemp))[0,1]
        zcrPerAxisGyro[j,2] = librosa.feature.zero_crossing_rate(GyroZtemp - statistics.mean(GyroZtemp))[0,1]
    
    # create dummy variables to fill in later with condition
    nameVectorPP = np.array([participants[p] for x in range(nEpochsOverlap)])

    # rewrite to panda dataframe format
    timeEpochs = pd.DataFrame(data=timeEpochs,columns=["Time"])
    nameVectorPP = pd.DataFrame(data=nameVectorPP,columns=["Participant"])
    Missingness = pd.DataFrame(data=np.array(["Valid" for x in range(nEpochsOverlap)]), columns=["Missingness"])
    # accelerometer
    averagePerAxis = pd.DataFrame(data=averagePerAxis,columns=["MeanX","MeanY","MeanZ","MeanVM"])
    stdPerAxis = pd.DataFrame(data=stdPerAxis,columns=["StdX","StdY","StdZ","StdVM"])
    variancePerAxis = pd.DataFrame(data=variancePerAxis,columns=["VarianceX","VarianceY","VarianceZ","VarianceVM"])
    maxPerAxis = pd.DataFrame(data=maxPerAxis,columns=["MaximumX","MaximumY","MaximumZ","MaximumVM"])
    maxNormPerAxis = pd.DataFrame(data=maxNormPerAxis,columns=["NormalizedMaximumX","NormalizedMaximumY","NormalizedMaximumZ","NormalizedMaximumVM"])
    minPerAxis = pd.DataFrame(data=minPerAxis,columns=["MinimumX","MinimumY","MinimumZ","MinimumVM"])
    minNormPerAxis = pd.DataFrame(data=minNormPerAxis,columns=["NormalizedMinimumX","NormalizedMinimumY","NormalizedMinimumZ","NormalizedMinimumVM"])
    rangePerAxis = pd.DataFrame(data=rangePerAxis,columns=["RangeX","RangeY","RangeZ","RangeVM"])
    kurtosisPerAxis = pd.DataFrame(data=kurtosisPerAxis,columns=["KurtosisX","KurtosisY","KurtosisZ","KurtosisVM"])
    skewnessPerAxis = pd.DataFrame(data=skewnessPerAxis,columns=["SkewnessX","SkewnessY","SkewnessZ","SkewnessVM"])
    corrPerAxis = pd.DataFrame(data=corrPerAxis,columns=["CorrXY","CorrXZ","CorrYZ"])
    inclPerAxis = pd.DataFrame(data=inclPerAxis,columns=["InclinationX","InclinationY","InclinationZ"])
    zcrPerAxis = pd.DataFrame(data=zcrPerAxis,columns=["zcrX","zcrY","zcrZ"])
    # gyroscope
    averagePerAxisGyro = pd.DataFrame(data=averagePerAxisGyro,columns=["MeanX_Gyro","MeanY_Gyro","MeanZ_Gyro","MeanVM_Gyro"])
    stdPerAxisGyro = pd.DataFrame(data=stdPerAxisGyro,columns=["StdX_Gyro","StdY_Gyro","StdZ_Gyro","StdVM_Gyro"])
    variancePerAxisGyro = pd.DataFrame(data=variancePerAxisGyro,columns=["VarianceX_Gyro","VarianceY_Gyro","VarianceZ_Gyro","VarianceVM_Gyro"])
    maxPerAxisGyro = pd.DataFrame(data=maxPerAxisGyro,columns=["MaximumX_Gyro","MaximumY_Gyro","MaximumZ_Gyro","MaximumVM_Gyro"])
    maxNormPerAxisGyro = pd.DataFrame(data=maxNormPerAxisGyro,columns=["NormalizedMaximumX_Gyro","NormalizedMaximumY_Gyro","NormalizedMaximumZ_Gyro","NormalizedMaximumVM_Gyro"])
    minPerAxisGyro = pd.DataFrame(data=minPerAxisGyro,columns=["MinimumX_Gyro","MinimumY_Gyro","MinimumZ_Gyro","MinimumVM_Gyro"])
    minNormPerAxisGyro = pd.DataFrame(data=minNormPerAxisGyro,columns=["NormalizedMinimumX_Gyro","NormalizedMinimumY_Gyro","NormalizedMinimumZ_Gyro","NormalizedMinimumVM_Gyro"])
    rangePerAxisGyro = pd.DataFrame(data=rangePerAxisGyro,columns=["RangeX_Gyro","RangeY_Gyro","RangeZ_Gyro","RangeVM_Gyro"])
    kurtosisPerAxisGyro = pd.DataFrame(data=kurtosisPerAxisGyro,columns=["KurtosisX_Gyro","KurtosisY_Gyro","KurtosisZ_Gyro","KurtosisVM_Gyro"])
    skewnessPerAxisGyro = pd.DataFrame(data=skewnessPerAxisGyro,columns=["SkewnessX_Gyro","SkewnessY_Gyro","SkewnessZ_Gyro","SkewnessVM_Gyro"])
    corrPerAxisGyro = pd.DataFrame(data=corrPerAxisGyro,columns=["CorrXY_Gyro","CorrXZ_Gyro","CorrYZ_Gyro"])
    inclPerAxisGyro = pd.DataFrame(data=inclPerAxisGyro,columns=["InclinationX_Gyro","InclinationY_Gyro","InclinationZ_Gyro"])
    zcrPerAxisGyro = pd.DataFrame(data=zcrPerAxisGyro,columns=["zcrX_Gyro","zcrY_Gyro","zcrZ_Gyro"])
    
    # concatenate into matrix
    featureMatrix = pd.concat([timeEpochs,nameVectorPP,Missingness,
                               averagePerAxis,stdPerAxis,variancePerAxis,maxPerAxis,maxNormPerAxis,minPerAxis,
                               minNormPerAxis,rangePerAxis,kurtosisPerAxis,skewnessPerAxis,corrPerAxis,inclPerAxis,zcrPerAxis,
                               averagePerAxisGyro,stdPerAxisGyro,variancePerAxisGyro,maxPerAxisGyro,maxNormPerAxisGyro,minPerAxisGyro,
                               minNormPerAxisGyro,rangePerAxisGyro,kurtosisPerAxisGyro,skewnessPerAxisGyro,corrPerAxisGyro,inclPerAxisGyro,zcrPerAxisGyro],axis=1)
    
    # rewrite time
    featureMatrix['Time'] = pd.to_datetime(featureMatrix['Time'], format="%Y-%m-%d %H:%M:%S", utc=True)
         
    if 'periods' in locals() or 'periods' in globals():
        for mp in range(len(periods)):
            startPeriod = pd.to_datetime(time[periods[mp][0]])
            endPeriod = pd.to_datetime(time[periods[mp][1]])
            millisecondsToStart = (timeEpochs - startPeriod).astype('timedelta64[ms]')
            millisecondsToEnd = (timeEpochs - endPeriod).astype('timedelta64[ms]')
            startMissing = np.where((millisecondsToStart > -3000) & (millisecondsToStart < 0))[0]
            endMissing = np.where((millisecondsToEnd > -3000) & (millisecondsToEnd < 0))[0]
            indexMissingEpoch = np.arange(startMissing, endMissing+1,dtype=int)
            # replace these indices with missing to avoid the epochs being used as training data
            featureMatrix.loc[indexMissingEpoch,'Missingness'] = 'Missing'
        del periods, startPeriod, endPeriod, millisecondsToStart, millisecondsToEnd, startMissing, endMissing, indexMissingEpoch
    
    # combine features into single matrix
    if p==0:
        postureFeatures = featureMatrix
    if p>0:
        postureFeatures = pd.concat([postureFeatures, featureMatrix], axis=0, ignore_index=True)

# save
postureFeatures.to_csv('Features.csv', index=False)
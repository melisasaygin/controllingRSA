# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 13:11:58 2024

@author: sve254
"""

# clear console and environment
#%reset -f 
import os
#os.system('cls')
import gc
#gc.collect()

# load libraries
import pandas as pd # handling dataframes
import numpy as np # array handling
import matplotlib.pyplot as plt
import seaborn as sns
import pickle

# %%# # create functions to transform epoch classification to window counts

def PredictionsToWindowCounts(predictions, TimeColumn, VMColumn, activities, epochsPerWindow):
    numWindows = int(len(predictions) // epochsPerWindow)
    countsPerWindow = pd.DataFrame(
        np.zeros((numWindows, len(activities)), dtype=int))
    countsPerWindow.columns = activities
    windowTime = []
    windowVM = []
    for w in range(numWindows):
        beginWindow = w*epochsPerWindow
        endWindow = beginWindow + epochsPerWindow
        windowWindow = predictions[beginWindow:endWindow]
        if len(windowWindow) == epochsPerWindow:
            windowTime.append(TimeColumn[beginWindow])
            windowVM.append(VMColumn[beginWindow:endWindow].mean())
            unique, counts = np.unique(windowWindow, return_counts=True)
            for u in range(len(unique)):
                countsPerWindow.loc[w, unique[u]] = counts[u]
    windowTime = pd.DataFrame(windowTime)
    windowTime.columns = ['Time']
    windowVM = pd.DataFrame(windowVM)
    windowVM.columns = ['MeanVM']
    windowTime = windowTime.join(windowVM)
    countsPerWindow = windowTime.join(countsPerWindow)
    return countsPerWindow

# NEW
def ProbabilitiesToWindowCounts(probabilities, TimeColumn, VMColumn, activities, epochsPerWindow):
    # find numer of windows and create empty
    numWindows = int(len(probabilities) // epochsPerWindow)
    probabilityPerWindow = pd.DataFrame(
        np.zeros((numWindows, len(activities)), dtype=int))
    probabilityPerWindow.columns = activities
    windowTime = []
    windowVM = []
    # loop over all windows and save counts and average probability
    for w in range(numWindows):
        # subset right part
        beginWindow = w*epochsPerWindow
        endWindow = beginWindow + epochsPerWindow
        probabilityWindow = probabilities[beginWindow:endWindow]
        if len(probabilityWindow) == epochsPerWindow:
            windowTime.append(TimeColumn[beginWindow])
            windowVM.append(VMColumn[beginWindow:endWindow].mean())
            # average probability per category
            averageProbability = np.mean(probabilityWindow, axis=0)
            probabilityPerWindow.iloc[w, :] = averageProbability
    # save
    windowTime = pd.DataFrame(windowTime)
    windowTime.columns = ['Time']
    windowVM = pd.DataFrame(windowVM)
    windowVM.columns = ['MeanVM']
    windowTime = windowTime.join(windowVM)
    probabilityPerWindow = windowTime.join(probabilityPerWindow)
    return probabilityPerWindow


def CombineCountsProbability(counts, probabilities):
    # create emtpy dataframe to save results
    majority = []
    for r in range(len(counts)):
        # subset row and find number of columns that contain max value
        row = counts.loc[r, :]
        max_value = row.max()
        max_columns = row[row == max_value].index
        if len(max_columns) == 1:
            majority.append(max_columns[0])
        elif len(max_columns) > 1:
           # get the probabilities for these rows
           probabilities_max_columns = probabilities.loc[r, max_columns]
           majority.append(probabilities_max_columns.idxmax())
    majority = np.array(majority)
    if len(majority) != len(counts):
        print(
            f'Different lengths counts and majority, please double check row {r}')
    return majority
# NEW

# set directory
os.chdir("C:/Users/msa583/Documents/New Posture Detections - Part 2/82556/day2")
folderpathData = "C:/Users/msa583/Documents/New Posture Detections - Part 2/82556/day2"

# read in data
Features = pd.read_csv('Features.csv', dtype={'Participant': str})

# find participants
participants = Features['Participant'].unique()

#%%# descriptive statistics full dataset

print(Features.shape)
print(Features.dtypes)
print(Features.describe())
print(Features.isnull().sum())
print(Features.iloc[:, -98:].corr())
    
# histograms
for col in Features.columns:
    if Features[col].dtypes == 'float64':
        sns.histplot(Features[col], bins=100)
        plt.title(col)
        plt.show()

#%%# load models and predict per participant

loadNameSVM = 'C:/Users/msa583/OneDrive - Vrije Universiteit Amsterdam/Desktop/Posture Detection Scripts/3. Feature extraction and prediction (in Spyder)/3. Feature extraction and prediction (in Spyder)/StaticDynamicClassifierSVM_Generic.pkl'
with open(loadNameSVM, 'rb') as f:
    StaticDynamicClassifierSVM_loaded = pickle.load(f)    

loadNameXGB = 'C:/Users/msa583/OneDrive - Vrije Universiteit Amsterdam/Desktop/Posture Detection Scripts/3. Feature extraction and prediction (in Spyder)/3. Feature extraction and prediction (in Spyder)/PostureClassifierXGB_Generic.pkl'
with open(loadNameXGB, 'rb') as f:
    PostureClassifierXGB_loaded = pickle.load(f)  
    
#%%# loop over all participants and predict category

participants = Features['Participant'].unique()

for p in range(len(participants)):
    print(participants[p])
    print(p+1, '/', len(participants))
    
    # subset data
    FeaturesTemp = Features.loc[Features['Participant'] == participants[p], :]

    #%%# static/dynamic classifier SVM
    
    # only use variables that are considered to be useful based on previous runs
    variableStrings = ['Time', 'Participant', 'Missingness', 'Std', 'Variance', 'Range', 'NormalizedMaximum']
    usefulColumns = [col for col in Features.columns if any(string in col for string in variableStrings)]
    
    # select useful columns
    FeaturesSVM = FeaturesTemp.loc[:, usefulColumns]
    FeaturesSVM = FeaturesSVM.iloc[:, 3:]
    
    # predict from model
    predictionSVM = StaticDynamicClassifierSVM_loaded.predict(FeaturesSVM)
    
    del variableStrings, usefulColumns
    
    #%%# posture detection XGBoost
    
    # filter out parts with unknown conditions
    FeaturesXGB = FeaturesTemp
    
    # drop irrelevant columns
    FeaturesXGB = FeaturesXGB.drop(columns=['Time', 'Participant', 'Missingness'])
    
    # predict from XGB
    predictionXGB = PostureClassifierXGB_loaded.predict(FeaturesXGB)
    probabilitiesXGB = PostureClassifierXGB_loaded.predict_proba(FeaturesXGB)
    
    # rewrite to names
    predictionXGB = pd.DataFrame(predictionXGB)
    predictionXGB.replace({0: 'Sitting', 1: 'Standing', 2: 'Lying'}, inplace=True)
    predictionXGB = predictionXGB.to_numpy().flatten()
    
    #%%# combine static/dynamic and postures
    
    combined_prediction = np.where(predictionSVM == 'Static', predictionXGB, predictionSVM)
    
    #%%# find majority within windows of certain length (30 seconds or 60 seconds)
    
    # set epoch length
    epochLength = 6
    overlap = 0
    windowLength = 60
    epochsPerWindow = int((windowLength/epochLength) * (1/(1-overlap)))
    
    # use function that was defined at the beginning of the script
    # HERE, CAN GET THE COUNTS PER POSTURE FOR EACH MINUTE! - save the window_counts_postures dataframe
    window_counts_postures = PredictionsToWindowCounts(predictionXGB, Features.Time, Features.MeanVM, ['Sitting', 'Standing', 'Lying'], epochsPerWindow)
    window_probabilities_postures = ProbabilitiesToWindowCounts(probabilitiesXGB, Features.Time, Features.MeanVM, ['Sitting', 'Standing', 'Lying'], epochsPerWindow)
    majority_postures = CombineCountsProbability(window_counts_postures[['Sitting', 'Standing', 'Lying']], window_probabilities_postures[['Sitting', 'Standing', 'Lying']])
    
    # use function that was defined 
    window_counts_including_dynamic = PredictionsToWindowCounts(combined_prediction, Features.Time, Features.MeanVM, ['Sitting', 'Standing', 'Lying', 'Dynamic'], epochsPerWindow)
    majority_including_dynamic = window_counts_including_dynamic[['Sitting', 'Standing', 'Lying', 'Dynamic']].idxmax(axis="columns").to_numpy()
    
    # create dataframe including time, participant and majorities
    majorities_temp = pd.concat([window_counts_including_dynamic['Time'], pd.DataFrame([participants[p]]*len(majority_postures), columns=['Participant']), pd.DataFrame(majority_postures, columns=['Majority_postures']), pd.DataFrame(majority_including_dynamic, columns=['Majority_including_dynamic'])], axis=1)
                       
    if p==0:
        majorities = majorities_temp
    if p>0:
        majorities = pd.concat([majorities, majorities_temp], axis=0)      

#%%# save majorities

majorities.to_csv('Majoriites per 60sec.csv', index=False)

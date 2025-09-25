# controllingRSA

This repository contains the code and data used in the study 'Controlling heart rate variability for respiratory effects in ambulatory psychophysiological measurements.' (doi upon publication). Please note that the dataset we provide here exclude the data of three participants who did not give consent for their pseudonymized data to be used in other studies. If you would like to replicate the results of the study to the exact statistics, and therefore need the full dataset, send an email to the corresponding author at m.saygin@vu.nl

There are 4 folders: Lab Calibration, Ambulatory Posture Detection, Merging Data and Calculating RSA Metrics, Analyses 

Lab Calibration:

Ambulatory Posture Detection: Includes: 1) the MATLAB script provided by the VU-AMS manufacturers (code authored by Aniket Mazumder and Cor Stoof of the VU Ambulatory Monitoring Solutions) to enable us extract the raw data, 2)Jupyter Notebook used to convert the .mat files to .h5, 3) the posture detection feature extraction and posture scripts as well as the machine learning model itself provided by Sjors van de Ven and colleagues.  

Merging Data and Calculating RSA Metrics:

Analyses: includes the R scripts used for multilevel model fitting and full variance decomposition of the models. The dataframes are also provided for all Research Questions 1, 2, and 3.

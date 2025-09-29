# controllingRSA

This repository contains the code and data used in the study 'Controlling heart rate variability for respiratory effects in ambulatory psychophysiological measurements.' (doi to be added upon publication). Please note that the dataset we provide here exclude the full data of three participants who did not give consent for their de-identified data to be used in other studies. To be absolutely sure we don't violate GDPR, we also removed the sensitive data (for all participants) including body mass index, affect/stress/safety self-reports, current smoking status, the physical activity level over the past week. If you are a researcher who would like the full dataset specifically to replicate the results of the study, please reac out to the corresponding author (m.saygin@vu.nl).

The de-identified physiological data of 39 participants are alreeady provided which should allow the analyses for Research Question 1 and 2, even without reaching out to the authors. While running the models for Research Question 2 in R, you would need to delete from code the covariates (BMI, smoking status, and PA are not included in the data) except for Age and Sex.

Organization of the Repository
There are 4 folders: Lab Calibration, Ambulatory Posture Detection, Merging Data and Calculating RSA Metrics, Analyses 

Lab Calibration:

Ambulatory Posture Detection: Includes: 1) the MATLAB script provided by the VU-AMS manufacturers (code authored by Aniket Mazumder and Cor Stoof of the VU Ambulatory Monitoring Solutions) to enable us extract the raw data, 2)Jupyter Notebook used to convert the .mat files to .h5, 3) the posture detection feature extraction and posture scripts as well as the machine learning model itself provided by Sjors van de Ven and colleagues.  

Merging Data and Calculating RSA Metrics:

Analyses: includes the R scripts used for multilevel model fitting and full variance decomposition of the models. The dataframes are also provided for all Research Questions 1, 2, and 3.

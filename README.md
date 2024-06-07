# Project CeLaVie
Scripts for the pre-processing and analysis of EEG data for the CeLaVie project.
### Data Preprocessing
The data preprocessing steps employ the **RELAX toolbox**, a Matlab-based toolbox that functions with **EEGLAB**.
Precisely, to detect ocular artifacts (vEOG, hEOG), muscle artifacts (EMG) and other noise sources such as drift, the RELAX toolbox employs a combination of **Multiple Wiener Filtering** and **Wavelet-enhanced Independent Components Analysis (ICA)**. 
The pipeline proposed by the RELAX toolbox also applies functions from the **PREP pipeline** (ref) to detect noisy electrodes and to apply **robust average referencing**. 

Below is presented an overview of the data preprocessing and analysis steps. 
<img title="Celavie data processing pipeline overview." alt="Alt text" src="/images/microstate_analysis_pipeline_flowchartv2.png">

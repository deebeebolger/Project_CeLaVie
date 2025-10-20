The function CLV_CreateRELAX() creates the RELAX_cfg.mat. This is the configuration structure containing the parameters that need to be defined prior to running the RELAX preprocessing pipeline. The following is the link to the CLV_CreateRELAX() function on Github.

### Important !
To be able to run the pipeline on your own computer, you will need to change the following three paths defined at the beginning of the CLV_CreateRELAX() function. 
1. **Relax_cfg.caploc** : Path to the channel location *.mat files
2. **Relax_cfg.myPath** : Path to the raw EEG data.
3. **Relax_cfg.OutputPath** : Path to the folder in which the processed data will be saved. This directory will be automatically created if not already present.

The configuration file is saved as a _*.json_ file for each subject and its title takes the following form :
 <span style="color:red;"> sub-XX_Relax_config.json </span>


## Downsampling : Calculating and recording parameters

Downsampling should, ideally, be carried out after low-pass filtering. However, given that the CeLaVie data was acquired at a sampling rate of 4096Hz and that a number of the preprocessing steps are computationally costly, it is necessary to reduce the number of samples before processing to other preprocessing steps. 

If the user decides to carry out downsampling, the parameter, ms par sample, needs to be calculated. This information is added to the RELAX configuration file.

Here the data is downsampled to a sampling rate of 512Hz and the ms per sample is calculated and added to the RELAX config file (RELAX_cfg) . 

```
RELAX_cfg.DownSample_to_X_Hz = 512
RELAX_cfg.ms_per_sample = 1000/RELAX_cfg.DownSample_to_X_Hz
```

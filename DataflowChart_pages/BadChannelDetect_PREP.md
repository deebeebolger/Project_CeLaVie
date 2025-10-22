## Noisy channel detection : PREP pipeline function

There are a total of three stages to the detection of noisy channels in the RELAX pipeline.
The automatic detection of the channels to be rejected, the following two criteria have to be defined and applied :
- In the GUI : *Maximum proportion of electrodes that can be deleted as bad* : In the parameters *RELAX_cfg.MaxProportionOfElectrodesThatCanBeDeleted* = 0.20 (20%)
- In the GUI : *Extreme noise proportion electrode deletion threshold* : RELAX_cfg.ProportionOfExtremeNoiseAboveWhichToRejectChannel = 0.05 (5%)

This first stage of noisy channel detection applies a function taken from the PREP pipeline (reference) : ```findNoisyChannels()```
This function divides the data into time windows of a pre-defined duration, a duration set to 1 second by default. Noisy channels are defined according to the following three criteria :

### Correlation with other channels - the *_predictability_* criterion
The signal is low pass filtered at 50Hz and divided into time windows of 1 second duration, and the correlation of each individual channel with other channels in each time window is calculated. The maximum absolute correlation is calculated as the 98th percentile of the absolute values of the correlations with other channels in each time window. A channel is marked as being bad by correlation if the maximum correlation is less than a defined correlation (default: 0.4) threshold for a given % of time windows (1%).  

### Unusual high frequency noise - the *_noisiness_* criterion
The ratio of the power of the high frequency components compared to the low frequency components is estimated. The ratio, for each channel, of the median absolute deviation (MAD) of high frequency components and  low frequency components is the noisiness. Based on this, the robust z-score is calculated for each channel compared to all other channels. Channels are marked as being bad-by-HF noise if their robust z-score > 5. This is calculated for each channel and for each time window. 

### Extreme amplitude - the *_deviation_* criterion
The deviation criterion is effective in detecting unusual high frequency noise. It is computed as the robust z-score of the robust standard deviation for each channel. Channels are marked as _bad-by-deviation_ if the robust z-score > 5. In addition, for each channel a robust amplitude-adjusted z score is calculated for each 1 second time window by using the overall robust median and the overall robust standard deviation in the z score calculation. 

## Filtering information

### Notch filtering
A 4th order Butterworth filter is applied with cutoff 3Hz above and below 50Hz, which implies low and high cutoff of 47Hz and 53Hz, respectively.


### Low-pass filtering
A 4th order Butterworth filter is applied with a high-pass cutoff of 0.25Hz. 
A lower cutoff ($\leq 0.1Hz$ ) would be preferable but given that we will be carrying out ICA and also given that we will not be carrying out ERP analysis on this data, we apply a cutoff of 0.25Hz; it has been shown that a higher high-pass cutoff can facilitate the computation of independent components.

### High-pass filtering
It is recommended not to carry out low pass filtering prior to Multi-channel Wiener Filtering (MWF) as this reduces the chance of rank deficiency. Therefore, in the current implementation of the pipeline, low-pass filtering is not carried at this stage. 

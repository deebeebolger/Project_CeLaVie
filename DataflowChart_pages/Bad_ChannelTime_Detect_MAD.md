## Noisy channel and time interval detection : MAD (Mean Absolute Deviation)

Extreme outlier detection is carried out using the *Median Absolute Deviation (MAD)*, which is more robust to outliers. The continuous data is segmented into 1 second epochs with a 50% overal (500ms) and, in each time segment, the following is identified *to
detect both electrodes and time segments* to mark for rejection :
- extreme amplitudes having a low probability of corresponding to brain activity. This is estimated by *MAD* from the median.
- extreme drift
- extreme kurtosis
- very improbable voltage distributions
- log-power log-frequency slopes, which can indicate muscle activity.

For each time window, a channel is deleted if the proportion of its data presenting extreme values exceeds a defined percentage, here **5%**. 
A maximum proportion of channels that can be deleted based on the above criterion is defined; here this proportion is 20%.


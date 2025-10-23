## Multi-channel Wiener filters : Round 2

In this second round of MWF, blinks are mainly dealt with. However, it is also possible to detect eye-blinks in the first round, in which case the eye-blink mask will be integrated into the noise mask. To allow eye-blink cleaning in round 2:
- **Relax_cfg.MWFRoundToCleanBlinks** = 2

As explained for the blink detection via IQR step, in the second round of MWF, the user needs to define the probability that the current data contains blink artifacts. The choices are:
- Data almost certainly has blinks : **Relax_cfg.ProbabilityDataHasNoBlinks** = 0
- Data might not have blinks : **Relax_cfg.ProbabilityDataHasNoBlinks** = 1
- Data definitely does not have blinks : **Relax_cfg.ProbabilityDataHasNoBlinks** = 2

In current implementation of the pipeline :
In GUI **_ Does the data contain blinks?_** : **Relax_cfg.ProbabilityDataHasNoBlinks** = 0

The aim is to detect eye-blinks that were masked my muscle activity and, thus, missed in round one of MWF.
To detect eye-blinks, a copy of data was bandpass filtered using a 4th order Butterworth filter between 1-25Hz. Data in channels identified as blink-relevant was averaged. And blinks were marked as the maximum point within each time-period that exceeded the value of the upper quartile of all voltages + the interquartile range (IQR)*3. 
This is carried out by the function in the following code if it is specified that the probability that the dataset does not have blinks is 0. 

```
[continuousEEG, epochedEEG] = RELAX_blinks_IQR_method(continuousEEG, epochedEEG,...
 RELAX_cfg)
```



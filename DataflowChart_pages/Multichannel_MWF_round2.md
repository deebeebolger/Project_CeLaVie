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
Thus, a blink mask was created for MWF cleaning by marking the 800ms surrounding all blink maximums as artifacts.

### Define blink-relevant channels

- In GUI **_Blink affected electrodes_** : **Relax_cfg.BlinkElectrodes** = C29, C17, C16, C30, C8, C28, C27, C21, C13, C10.
- In GUI **_Left sided HEOG affected electrodes_** : **Relax_cfg.HEOGLeftpattern** = C30, D7, C8, D9, D23, D10, D22, D24, C31.
- In GUI **_Right-sided HEOG affected electrodes_** : **Relax_cfg.HEOGRightpattern** = C8, C7, B27, B28, B26, B29, B24, B14, C9.

 ![Eyeblink channels](https://github.com/deebeebolger/Project_CeLaVie/blob/main/DataflowChart_pages/BlinkChannels_figure.png)
 Topopgraphy showing the locations of the defined blink channels passed to the blink detection that is carried out using the IQR routine.



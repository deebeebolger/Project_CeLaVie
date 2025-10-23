## Multi-channel Wiener filter : Round 3

In this 3rd round of Multi-channel Wiener Filtering (MWF), drift and horizontal ocular artifacts are detected.  The user has the option of running this 3rd MWF round.

To detect drift, the continuous EEG is re-referenced using PREP’s robust average referencing approach and epochs showing an amplitude at any electrode greater than a threshold were marked as artifact periods in the template used for MWF cleaning. Note that the re-referencing is only applied to identify drift and the re-referenced data is not used in the MWF cleaning. 

Drift detection is carried out with the following function :
```
[continuousEEG, epochedEEG] = RELAX_drift(continuousEEG, epochedEEG, RELAX_cfg)
```
The drift threshold is the mean absolute deviation (MAD) from the median of all electrodes. 
- In GUI **_Single electrode drift threshold for MWF cleaning_°° : **Relex_cfg.DriftSeverityThreshold** = 10

The maximum proportion of epochs to include in the mask from the drift artifact type.
- In GUI **_Max proportion marked as drift for MWF cleaning_°° : **Relex_cfg.ProportionWorstEpochsForDrift** = 0.3

The horizontal eye-movement threshold calculated as the MAD deviation from the median. if both lateral electrodes show activity above this threshold for a certain duration ( the number of timepoints (ms) that exceed the horizontal eye-movement threshold), the activity is marked as horizontal eye-movement. 

- In GUI **_Horizontal eye movement threshold (MAD from median) for MWF cleaning_** : **Relex_cfg.HorizontalEyeMovementThreshold** = 2




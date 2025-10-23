## IC classification and wavelet thresholding of ICs

The computed ICs (Independent Components) are automatically classified using EEGLAB’s ICLabel plugin (reference). The algorithms underlying the automatic classification of ICs is based on crowd sourcing run by the Swartz Center for Computational Neuroscience, University of California, which involves experts contributing to the ICLabel project by volunteering to label components. The ICLabel plugin automatically classifies ICs into broad categories based on their source:
- Brain
- Muscle
- Ocular
- Cardiac
- Channel noise
- Line noise
- Other

In the current implementation, an IC is considered for rejection if the probability of it corresponding to **Muscle**, **Ocular**, **Cardiac**, **Channel or Line** noise is equal to or above 80%. 

```
% IC classes: Brain, Muscle, Eye, Heart, LineNoise, ChannelNoise, Other
 icThreshold    = [0 0.0; 0.8 1; 0.8 1; 0.8 1; 0.8 1; 0.8 1; 0 0];
```
Matlab code defining the criteria for marking for rejection those ICs classified by ICLabel.






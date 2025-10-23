## Multi-channel Wiener filters : Round 2

In this second round of MWF, blinks are mainly dealt with. However, it is also possible to detect eye-blinks in the first round, in which case the eye-blink mask will be integrated into the noise mask. To allow eye-blink cleaning in round 2:
- **Relax_cfg.MWFRoundToCleanBlinks** = 2

As explained for the blink detection via IQR step, in the second round of MWF, the user needs to define the probability that the current data contains blink artifacts. The choices are:
- Data almost certainly has blinks : **Relax_cfg.ProbabilityDataHasNoBlinks** = 0
- Data might not have blinks : **Relax_cfg.ProbabilityDataHasNoBlinks** = 1
- Data definitely does not have blinks : **Relax_cfg.ProbabilityDataHasNoBlinks** = 2

In current implementation of the pipeline :
In GUI **_ Does the data contain blinks?_** : **Relax_cfg.ProbabilityDataHasNoBlinks** = 0


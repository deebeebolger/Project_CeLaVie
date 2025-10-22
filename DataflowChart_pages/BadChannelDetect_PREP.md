## Noisy channel detection : PREP pipeline function

There are a total of three stages to the detection of noisy channels in the RELAX pipeline.
The automatic detection of the channels to be rejected, the following two criteria have to be defined and applied :
- In the GUI : *Maximum proportion of electrodes that can be deleted as bad* : In the parameters *RELAX_cfg.MaxProportionOfElectrodesThatCanBeDeleted* = 0.20 (20%)
- In the GUI : *Extreme noise proportion electrode deletion threshold* : RELAX_cfg.ProportionOfExtremeNoiseAboveWhichToRejectChannel = 0.05 (5%)

This first stage applies a function taken from the PREP pipeline (reference) : ```*_findNoisyChannels()_*```

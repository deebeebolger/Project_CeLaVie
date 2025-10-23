## Multi-channel Wiener filter : Round 3

In this 3rd round of Multi-channel Wiener Filtering (MWF), drift and horizontal ocular artifacts are detected.  The user has the option of running this 3rd MWF round.

To detect drift, the continuous EEG is re-referenced using PREP’s robust average referencing approach and epochs showing an amplitude at any electrode greater than a threshold were marked as artifact periods in the template used for MWF cleaning. Note that the re-referencing is only applied to identify drift and the re-referenced data is not used in the MWF cleaning. 

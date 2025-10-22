## Eye-blink detection : IQR approach

This method creates a mask that specifies the localisation of eye-blinks in continuous EEG by using the inter-quartile range threshold.

>| :memo:        | This method is applied in the second round of Multi-channel Wiener Filtering (MWF) and so, if at least the first and second rounds of MWF are being carried out, it is not necessary to carry out this step here.       |

In the **RELAX config**, the user specifies the probability that the data does not have blinks, as follows :
- 0 $\implies$ data almost certainly has blinks
- 1 $\implies$ data might not have blinks
- 2 $\implies$ data definitely does not have blinks






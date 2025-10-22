## Eye-blink detection : IQR approach

This method creates a mask that specifies the localisation of eye-blinks in continuous EEG by using the inter-quartile range threshold.
[!NOTE] 
This method is applied in the second round of Multi-channel Wiener Filtering (MWF) and so, if at least the first and second rounds of MWF are being carried out, it is not necessary to carry out this step here.


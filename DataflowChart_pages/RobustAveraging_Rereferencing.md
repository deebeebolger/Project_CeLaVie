## Robust average re-referencing

After Multi-channel Wiener Filter (MWF) stages, the continuous data is re-referenced to the average across the scalp channels, excluding those channels marked as bad. 

|❗: In the RELAX re-reference function, RELAX_average_rereference(EEG), is described as applying the robust re-reference but this does not appear to be the case. This function applies EEGLAB's ```pop_reref()``` function, which re-references to the average of the scalp electrodes.To apply the PREP pipeline's robust re-referencing function, it is necessary to use their ```performReference()```or ```robustReference()``` function. |
|--------------------------------------------------------------------------------------------------------------|

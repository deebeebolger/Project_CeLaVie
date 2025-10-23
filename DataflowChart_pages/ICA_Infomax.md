## ICA via le Infomax algorithm

Independent Components Analysis (ICA) is a statistical and computational technique for blind source separation and feature extraction or to express this otherwise, for uncovering the hidden factors underlying sets of signals, random variables or measurements. The data variables that we feed into the ICA are assumed to be linear mixtures of a set of **_latent variables_** and both the latent variables and the mixing system is unknown. Importantly, the latent variables that ICA reveals, which we generally refer to as **_independent components_**, are assumed to be both **_non-gaussian_** and **_mutually independent_**. ICA is related to PCA, with the main difference that PCA assumes orthogonality and not independence, thus rendering it less powerful. PCA is often used in conjunction with ICA to carry out dimension reduction prior to ICA computation. 

Independent Components Analysis (ICA) is carried out after Multi-channel Wiener Filtering to detect and correct noise that MWF has difficulty dealing with. MWF has the disadvantage of requiring templates of both artifacts and clean data, which presents a difficulty when dealing with artifacts that occurs consistently throughout the data such as persistent muscle artifacts or electrical interference (line noise). ICA, on the other, is used to deal with the problem of blind source separation; when we do not have an example of the noise sources that we are trying to separate. ICA is particularly suited to neural signals as it does not require the assumption of normality, which is case with PCA (principle components analysis). 

### Infomax Algorithm

To apply ICA, the InfoMax (information maximisation) algorithm (Bell & Sejnowski, 1995a, 1995b) is applied using EEGLAB’s ```infomax()`` function.  Specifically, we apply the extended version of the InfoMax algorithm, which can deal with noise super-Gaussian and sub-Gaussian noise sources. This algorithm is more computationally expensive but gives more reliable results. Channels marked as bad are not included in the ICA and the number of ICs (independent components) computed equals the number of channels (excluding bad channels marked for rejection). 

Artifacts detected by ICA are considered as being independent of neural activity (our signal of interest) and can theoretically, therefore, be subtracted from the signal. However, ICA is very sensitive to high amplitude noise sources; their presence adversely affects the ICA decomposition leading to ICs that are a mixture of both neural and noise sources. This implies the danger of eliminating both noise and neural activity when carrying out noise correction with ICA.  It is for this reason that noisy time periods and bad channels should not be included in the ICA decomposition. However, even in the absence of very large amplitude artifacts, ICA can return components that are a mixture of neural and noise signals, although generally such components account for a smaller proportion of the data than those ICs that have successfully captured artifacts only. 

it is also for this reason that **wavelet-enhanced ICA (wICA)** is carried out in this pipeline.

```
% iseeg ; variable with the indices of good scalp electrodes.
% ncomps ; the number of ICs, which is equal to the number of good, scalp channels.

[weights, sphere] = runica(EEG.data(iseeg,:), 'ncomps', length(ncomps), 'extended', 1); 

%% Copy the IC weigths and sphere information to EEG dataset.

EEG.icaweights = weights;
EEG.icasphere  = sphere;
EEG.icachansind = iseeg;
W = weights * sphere;
A = inv(W);
EEG.icawinv = A;
EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(iseeg,:); % Calculate the ICA activations.
```
Matlab code in pipeline to run ICA (extended InfoMax algorithm) and calculate the ICs.


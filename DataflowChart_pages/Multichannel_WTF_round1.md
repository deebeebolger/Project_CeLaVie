## Multi-channel Wiener filters : Round 1

Multi-channel Wiener Filters (MWF) (Borowicz, 2018; Somers et al., 2018) constitutes one of the key artifact reduction components of RELAX pipeline. The MWF are applied sequentially in three rounds to reduce:
1. Muscle artifacts
2. Blinks that may have been masked by muscle artifacts
3. Horizontal eye-movements (hEOG) and drift

The aim, in this step, is to reduce as much as possible these artifact-types while preserving, as much as possible, the neural signal.

### Background
Wiener filters are data-dependent, linear, least square error filters. The coefficients of Wiener filters are calculated so as to minimize the average square distance between the filter output and a desired signal. In its most basic form, Wiener filters assume that signals are stationary processes but by periodically recalculating the filter coefficients for every block of N signal samples, the filter can adapt itself to the average characteristics of the signals within the block. 

Here, the activity underlying the EEG signal cannot be considered stationary and so blocks of N samples need to be defined. In the pipeline, these N samples define the **MWF delay period** and need this delay period is calculated as a function of the sampling frequency. 

If the desired delay period (D) is 40ms (0.04seconds) then N (the MWF delay period in samples) is given by $N = floor(D * F_s)$ where $F_s$ is the sampling frequency (512Hz).

In the current implementation of the pipeline, a delay period of 40 milliseconds is defined, which translates as a delay period (in samples) of **20 samples**.

In GUI **_MWF Delay Period_** : **Relax_cfg.MWFDelayPeriod_ms** = 40

| :exclamation:  For a delay period of 16 samples, we would need to defined a delay period (in ms) of 31.2ms   |
|--------------------------------------------------------------------------------------------------------------|

### Multi-channel Wiener Filtering (MWF) : Round 1

The **first round** concerns mainly **_muscle artifact_** cleaning. But eye-blink cleaning can be included in this first round if wished. This implies integrating the eye-blink mask into the noise mask. The first round is carried if : **Relax_cfg.Do_MWF_Once** = 1

#### Muscle activity

Muscle artifacts are dealt with in the first round of MWF. To create a template for the detection of muscle activity for the MWF cleaning, the data is separated into 1second epochs with 500ms overlap.



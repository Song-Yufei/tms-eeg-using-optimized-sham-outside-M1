# tms-eeg-using-optimized-sham-outside-M1
This repository contains the code used to produce results in the manuscript: Song, Y., Gordon, P. C., Metsomaa, J., Rostami, M., Belardinelli, P., & Ziemann, U. (2023). Evoked EEG Responses to TMS Targeting Regions Outside the Primary Motor Cortex and Their Test–Retest Reliability. Brain Topography, 1-18. 

Summary:

We aimed to investigate TEPs and their test-retest reliability when targeting regions outside M1, specifically the left angular gyrus (AG), supplementary motor area (SMA), and medial prefrontal cortex (mPFC), using an optimized sham procedure. We conducted three identical TMS–EEG sessions one week apart involving 24 healthy participants. In each session, we targeted the three areas separately using a figure−of−eight TMS coil for active TMS, while a second coil away from the head produced auditory input for sham TMS. Masking noise and electric scalp stimulation were applied in both conditions to achieve matched EEG responses to peripheral sensory inputs. 
TMS-EEG recording:

For each TMS-EEG session, we divided it into three TMS-EEG blocks that correspond to three cortical targets: AG, SMA, and mPFC. EEG signals were recorded with a TMS-compatible system (NeurOne, Bittium). Electrodes were placed according to the International 10-5 system in an elastic cap (EasyCap BC-TMS-64, EasyCap). EEG was sampled at 5 kHz (device filter DC-1250Hz), and electrode CPz served as the reference online. We recorded 150 pulses for the active conditions per cortical target and 150 for sham conditions.

Visual analog scale (VAS):

Participants rated their perception of auditory and somatosensory inputs from the active and sham conditions, ranging from 0 to 10. 0 represented no perception, and 10 described maximal perception. The VAS included items assessing the intensity of auditory sensation, the intensity of scalp sensation, the area size of scalp sensation, and the intensity of pain or discomfort. Each item rating was replicated twice.

Dataset:

The preprocessed eeg data (in fieldtrip data structure), VAS scores and session information are deposited at the repository Zenodo.

Data analysis:

Offline EEG analyses were performed in Matlab environment. EEGLAB, FieldTrip, RStudio and customized scripts were used. 

‘tms_eeg_Cleaningpipeline.m’ lists the preprocessing steps for TMS-EEG data. 

‘load_grandAvg_tep.m’ and ‘load_grandAvg_gmfa.m’ list the steps: loads the preprocessed data after the cleaning pipeline, calculates stimulus-related averaging and GMFA per participant per cortical target, and makes a grand average over participants.

‘spatialCCC.m and temporalCCC.m lists the steps for intersession concordance correlation coefficient (CCC) (S1 v S2, S1 v S3, S2 v S3) in spatial and temporal domains per cortical target. Temporal CCCs were calculated within each TOI and the whole-time window. Peak latency and range were determined using ‘findpeaks.m’ and saved as ‘peak_latency_ranges_ag.mat’, ‘peak_latency_ranges_sma.mat’ and ‘peak_latency_ranges_mpfc.mat’. One-sample permutation t-tests were used for statistical analysis at the group level. ‘individualCCC.m’ calculates and illustrates CCCs at the individual level. 

‘statistics_tep.m’ lists the steps for statistically comparing the EEG spatiotemporal profile between the active and sham conditions per cortical target at the group level. The cluster-based permutation t-test was applied to compare evoked EEG potentials from both conditions across time points and electrodes within each TOI. 

‘reconstructSource_tep.m’ lists the steps: loads the preprocessed data after the cleaning pipeline, individual headmodel and leadfield matrix, and estimates source activity of TEPs per participant per cortical target. The reconstructed signal was z-score normalized with respect to a pre-stimuli time window (-600 to -100 ms).

‘VAS_sensation.R’ loads VAS scores and performs the Wilcoxon signed-rank test.

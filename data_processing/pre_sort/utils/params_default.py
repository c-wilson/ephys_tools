"""Contain default parameters."""

assert sample_rate > 0

# Filtering
# ---------
filter_low = 500.  # Low pass frequency (Hz)
filter_high = 0.95 * .5 * sample_rate
filter_butter_order = 3  # Order of Butterworth filter.
pl_trig_chan = None  # SpikeGL channel number where the PL trigger resides.

# Chunks
# ------
chunk_size = int(1. * sample_rate)  # 1 second
chunk_overlap = int(.015 * sample_rate)  # 15 ms

# Saving raw/filtered edata
# ------------------------
save_raw = True
save_high = False
save_low = True

# Spike detection
# ---------------
# Uniformly scattered chunks, for computing the threshold from the std of the
# signal across the whole recording.
nexcerpts = 50
excerpt_size = int(1. * sample_rate)
threshold_strong_std_factor = 4.5
threshold_weak_std_factor = 2.
detect_spikes = 'negative'

# Connected component
# -------------------
connected_component_join_size = int(.00005 * sample_rate)

# Spike extraction
# ----------------
extract_s_before = 10
extract_s_after = 10
waveforms_nsamples = extract_s_before + extract_s_after

# Features
# --------
nfeatures_per_channel = 3  # Number of features per channel.
pca_nwaveforms_max = 10000
features_contiguous = True  # Whether to make the features array contiguous



# ########################
# KlustaKwik parameters #
#########################
MaskStarts = 100
#MinClusters = 100 
#MaxClusters = 110
MaxPossibleClusters = 500
FullStepEvery = 10
MaxIter = 10000
RandomSeed = 654
Debug = 0
SplitFirst = 20
SplitEvery = 100
PenaltyK = 0
PenaltyKLogN = 1
Subset = 1
PriorPoint = 1
SaveSorted = 0
SaveCovarianceMeans = 0
UseMaskedInitialConditions = 1
AssignToFirstClosestMask = 1
UseDistributional = 1



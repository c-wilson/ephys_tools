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
kk_MaskStarts = 500
#kk_MinClusters = 100
#kk_MaxClusters = 110
kk_MaxPossibleClusters = 500
kk_FullStepEvery = 10
kk_MaxIter = 10000
kk_RandomSeed = 654
kk_Debug = 0
kk_SplitFirst = 20
kk_SplitEvery = 40
kk_PenaltyK = 1
kk_PenaltyKLogN = 0
kk_Subset = 1
kk_PriorPoint = 1
kk_SaveSorted = 0
kk_SaveCovarianceMeans = 0
kk_UseMaskedInitialConditions = 1
kk_AssignToFirstClosestMask = 1
kk_UseDistributional = 1
kk_RamLimitGB = 120.


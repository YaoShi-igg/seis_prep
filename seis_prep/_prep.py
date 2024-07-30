
def prep(trace, dest_sampling_rate, freq_range=None):
    """
    Prepare the seismic trace for analysis by resampling, filtering, and detrending.
    Args:
        - trace (obspy.core.trace.Trace): The input seismic trace to be processed.
        - dest_sampling_rate (float): The desired sampling rate after resampling.
        - freq_range (tuple, optional): A tuple of (freqmin, freqmax) for bandpass filtering.
    Returns:
        - obspy.core.trace.Trace: The processed seismic trace.
    """
    tr = trace.copy()
    if tr.stats.sampling_rate != dest_sampling_rate:
        tr.resample(sampling_rate=dest_sampling_rate, window='hann', no_filter=True, strict_length=False)
    
    # must before filter !!!!
    tr.detrend('demean')
    tr.detrend('linear')

    if freq_range is not None:
        tr.filter("bandpass", freqmin=freq_range[0], freqmax=freq_range[1])

    return tr

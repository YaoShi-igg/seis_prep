import numpy as np
import math
from matplotlib.colors import Normalize
from matplotlib import mlab
from obspy.imaging.cm import obspy_sequential

def _nearest_pow_2(x):
    """
    Find the power of two nearest to x.
    Args:
        - x (float): Number to find the nearest power of two for.
    Returns:
        - int: Nearest power of 2 to x.
    """
    if x < 1:
        return 1
    
    upper_pow = 2 ** math.ceil(math.log2(x))
    lower_pow = 2 ** math.floor(math.log2(x))
    
    if abs(upper_pow - x) < abs(lower_pow - x):
        return upper_pow
    else:
        return lower_pow

def spectrogram(data, samp_rate, axes, per_lap=0.9, wlen=None, 
                dbscale=False, mult=8.0, cmap='viridis', 
                zorder=None, clip=[0.0, 1.0]):
    """
    Computes and plots the spectrogram of the input data.
    Args:
        - data (np.ndarray): Input data array.
        - samp_rate (float): Sampling rate in Hz.
        - axes (matplotlib.axes.Axes): Axes object to plot the spectrogram.
        - per_lap (float): Percentage of overlap of sliding window, ranging from 0 to 1.
        - wlen (int or float): Window length for FFT in seconds. Defaults to 128 samples if None.
        - dbscale (bool): If True, use 10 * log10 of color values, else use sqrt.
        - mult (float): Pad zeros to length mult * wlen to make the spectrogram smoother.
        - cmap (str or matplotlib.colors.Colormap): Colormap for the spectrogram.
        - zorder (float): Z-order for the plot.
        - clip (list of float): Adjust colormap to clip at lower and/or upper end.
    Returns:
        - None
    """
    # Enforce float for samp_rate
    samp_rate = float(samp_rate)

    # Set default window length if not specified
    if wlen is None:
        wlen = 128 / samp_rate
    npts = len(data)

    # Calculate the FFT length
    nfft = int(_nearest_pow_2(wlen * samp_rate))

    if npts < nfft:
        raise ValueError(f'Input signal too short ({npts} samples) for the given window length ({wlen} seconds, nfft {nfft} samples) and sampling rate ({samp_rate} Hz).')

    if mult is not None:
        mult = int(_nearest_pow_2(mult))
        mult = mult * nfft
    nlap = int(nfft * float(per_lap))

    # Remove mean to avoid DC offset
    data = data - data.mean()

    # Compute the spectrogram
    specgram, freq, time = mlab.specgram(data, Fs=samp_rate, NFFT=nfft, pad_to=mult, noverlap=nlap)

    if len(time) < 2:
        raise ValueError(f'Input signal too short ({npts} samples) for the given window length ({wlen} seconds, nfft {nfft} samples) and window overlap ({nlap} samples) with sampling rate ({samp_rate} Hz).')

    # Apply db scale or sqrt scale
    if dbscale:
        specgram = 10 * np.log10(specgram[1:, :])
    else:
        specgram = np.sqrt(specgram[1:, :])
    freq = freq[1:]

    # Apply clipping to the spectrogram
    vmin, vmax = clip
    if vmin < 0 or vmax > 1 or vmin >= vmax:
        raise ValueError("Invalid parameters for clip option.")
    
    _range = float(specgram.max() - specgram.min())
    vmin = specgram.min() + vmin * _range
    vmax = specgram.min() + vmax * _range
    norm = Normalize(vmin, vmax, clip=True)

    # Calculate half bin width for centering
    halfbin_time = (time[1] - time[0]) / 2.0
    halfbin_freq = (freq[1] - freq[0]) / 2.0

    # Prepare arguments for pcolormesh
    kwargs = {k: v for k, v in (('cmap', cmap), ('zorder', zorder)) if v is not None}

    # Adjust the frequency and time arrays for pcolormesh
    freq = np.concatenate((freq, [freq[-1] + 2 * halfbin_freq]))
    time = np.concatenate((time, [time[-1] + 2 * halfbin_time]))
    time -= halfbin_time
    freq -= halfbin_freq

    # Plot the spectrogram
    cax = axes.pcolormesh(time, freq, specgram, norm=norm, rasterized=True, **kwargs)
    axes.set_yscale('log')
    axes.grid(False)

    # plt.colorbar(cax, ax=axes, label='Amplitude')

    return

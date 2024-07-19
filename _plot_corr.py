import os
import glob
import numpy as np
import math
from obspy import Stream
from obspy import read
import matplotlib.pyplot as plt
try:
    from ._cal_func import _nearest_pow_2, spectrogram
except:
    from _cal_func import _nearest_pow_2, spectrogram

def plot_single_ncf(ncf_path, filter_range, stft, fig_dir, fig_type):
    """
    Plots a single noise correlation function (NCF) with its corresponding spectrogram.
    Args:
        - ncf_path (str): Path to the NCF file.
        - filter_range (tuple): Frequency range for bandpass filter (min_freq, max_freq). Use None for no filter.
        - stft (Bool): If plot the Time-frequency figures for ncf.
        - fig_dir (str): Directory to save the figure.
        - fig_type (str): File type for the figure (e.g., 'png', 'jpg').
    Returns:
        - None
    """
    # Read the NCF file
    tr = read(ncf_path)[0]
    stack_type = ncf_path.split('.')[-1]
    title = ('.'.join(ncf_path.split('.')[:-1])).split('/')[-1]

    # Apply bandpass filter if specified
    if filter_range is not None:
        tr.filter('bandpass', freqmin=filter_range[0], freqmax=filter_range[1])
        title = f"{title} - {round(tr.stats.sac.dist)} km - {filter_range[0]}-{filter_range[1]}Hz"
    else:
        title = f"{title} - {round(tr.stats.sac.dist)} km"

    # Get data and time information
    len_ncf = tr.stats.endtime - tr.stats.starttime
    dt = tr.stats.delta
    t = np.arange(-len_ncf/2, len_ncf/2 + dt, dt)
    data = tr.data

    if stft:
        # Create the figure and axes
        fig, ax = plt.subplots(2, 1, figsize=(10, 6))
        # Plot the NCF
        ax[0].plot(t, tr.data, color='black', linewidth=0.5)
        ax[0].axvline(x=0, color='red')
        ax[0].set_xlabel('Lag Time (s)', fontsize=15)
        ax[0].set_ylabel('Normalized Amplitude', fontsize=15)
        ax[0].set_title(title, fontsize=16)
        ax[0].set_xlim(-len_ncf/2, len_ncf/2)
        ax[0].set_xticks(np.linspace(-len_ncf/2, len_ncf/2, 9))
        ax[0].set_xticklabels(np.linspace(-len_ncf/2, len_ncf/2, 9).astype(int))
        ax[0].xaxis.set_tick_params(labelsize=13)
        ax[0].yaxis.set_tick_params(labelsize=13)

        # Plot the spectrogram
        spectrogram(data=data, samp_rate=tr.stats.sampling_rate, axes=ax[1], cmap='hot')
        ax[1].set_xlim(0, len_ncf)
        ax[1].set_ylabel('Frequency (Hz)', fontsize=15)
        ax[1].set_xticks(np.linspace(0, len_ncf, 9))
        ax[1].set_xticklabels(np.linspace(0, len_ncf, 9).astype(int))
        ax[1].xaxis.set_tick_params(labelsize=13)
        ax[1].yaxis.set_tick_params(labelsize=13)

        fig_name = f"{title}.stft.{fig_type}"
        plt.tight_layout()

    else:
        fig, ax = plt.subplots(1, 1, figsize=(10, 5))
        # Plot the NCF
        ax.plot(t, tr.data, color='black', linewidth=0.5)
        ax.axvline(x=0, color='red')
        ax.set_xlabel('Lag Time (s)', fontsize=15)
        ax.set_ylabel('Normalized Amplitude', fontsize=15)
        ax.set_title(title, fontsize=16)
        ax.set_xlim(-len_ncf/2, len_ncf/2)
        ax.set_xticks(np.linspace(-len_ncf/2, len_ncf/2, 9))
        ax.set_xticklabels(np.linspace(-len_ncf/2, len_ncf/2, 9).astype(int))
        ax.xaxis.set_tick_params(labelsize=13)
        ax.yaxis.set_tick_params(labelsize=13)
        fig_name = f"{title}.{fig_type}"

    fig_path = os.path.join(fig_dir, fig_name)
    plt.savefig(fig_path, dpi=300)
    plt.close()

    return

def plot_ncfs(ncf_path_list, gain, filter_range, fig_dir, fig_type):
    """
    Plots a list of noise correlation functions (NCFs).
    Args:
        - ncf_path_list: List of paths to NCF files.
        - gain: Gain to scale the NCF data.
        - filter_range: Tuple (freqmin, freqmax) specifying the bandpass filter range.
        - fig_dir: Directory to save the generated plots.
        - fig_type: File format of the saved plots (e.g., 'png', 'pdf').
    Returns:
        - None
    """
    # Extract station name from the first file path
    center_sta = (ncf_path_list[0].split('/')[-1]).split('_')[0]
    # print(center_sta)

    # Read all NCF data into a single Stream object
    st = Stream()
    for ncf_path in ncf_path_list:
        st += read(ncf_path)

    # Apply bandpass filter if specified
    if filter_range is not None:
        st.filter('bandpass', freqmin=filter_range[0], freqmax=filter_range[1])    

    # params for plot
    len_ncf = st[0].stats.endtime - st[0].stats.starttime
    dt = st[0].stats.delta
    t = np.arange(-len_ncf/2, len_ncf/2 + dt, dt)

    # plot the ncfs
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    fig_path = os.path.join(fig_dir, '.'.join([center_sta, 'ncfs', fig_type]))

    ax.set_xlabel("Lag Time (s)", fontsize=20)
    ax.set_ylabel("Distance (km)", fontsize=20)
    ax.set_xlim(-len_ncf/2, len_ncf/2)
    ax.xaxis.set_tick_params(labelsize=18)
    ax.yaxis.set_tick_params(labelsize=18)

    if filter_range is not None:
        ax.set_title(f"{center_sta}_{filter_range[0]}-{filter_range[1]}Hz", fontsize=22)
    else:
        ax.set_title(f"{center_sta}_nofilter", fontsize=22)

    for tr in st:
        data = tr.data / max(tr.data) * gain + tr.stats.sac['dist']
        ax.plot(t, data, color='black', linewidth=0.5)

    plt.savefig(fig_path, dpi=300)
    plt.close(fig)

    return

if __name__ == '__main__':

    ncf_path = './ncf_data/AEJ.1001_AEJ.1005_ls.SAC'
    ncf_path_list = glob.glob('./ncf_data/*')
    # print(ncf_path_list)

    plot_single_ncf(ncf_path=ncf_path, filter_range=[0.1, 1], stft=True, fig_dir='./figures', fig_type='png')
    plot_single_ncf(ncf_path=ncf_path, filter_range=None, stft=False, fig_dir='./figures', fig_type='png')

    plot_ncfs(ncf_path_list=ncf_path_list, gain=3, filter_range=None, fig_dir='./figures', fig_type='png')

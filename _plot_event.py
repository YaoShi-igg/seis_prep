import os
import glob
from obspy import UTCDateTime
import numpy as np
import math
from obspy import read, Stream
import matplotlib.pyplot as plt
try:
    from ._time import utc2int
    from ._cal_func import _nearest_pow_2, spectrogram
except:
    from _time import utc2int
    from _cal_func import _nearest_pow_2, spectrogram

def plot_3chns_event(filepath_chn1, filepath_chn2, filepath_chn3, filter_range, figure_dir, figure_type):
    """
    Plot seismograms for 3 channels.
    Args:
        - filepath_chn1 (str): Path to the first channel data file.
        - filepath_chn2 (str): Path to the second channel data file.
        - filepath_chn3 (str): Path to the third channel data file.
        - filter_range (tuple): Tuple containing the minimum and maximum frequencies for bandpass filter.
        - figure_dir (str): Directory to save the figure.
        - figure_type (str): Type of the figure file (e.g., 'png', 'jpg').
    Returns:
        - None
    """
    plt.rcParams['xtick.labelsize'] = 12  # Set the font size for tick labels on the x-axis
    plt.rcParams['ytick.labelsize'] = 12
    # read the traces
    tr1, tr2, tr3 = read(filepath_chn1)[0], read(filepath_chn2)[0], read(filepath_chn3)[0]

    # bandpass filter
    if filter_range != None:

        tr1.filter('bandpass', freqmin=filter_range[0], freqmax=filter_range[1])
        tr2.filter('bandpass', freqmin=filter_range[0], freqmax=filter_range[1])
        tr3.filter('bandpass', freqmin=filter_range[0], freqmax=filter_range[1])
    
    # read the data
    data1, data2, data3 = tr1.data, tr2.data, tr3.data
    sampling_rate = tr1.stats.sampling_rate
    t = np.arange(0, len(data1)) / sampling_rate

    # P arrival location on the waveform
    p_loc = tr1.stats.sac.t1 - tr1.stats.sac.b

    # Plot the figure
    fig, ax = plt.subplots(3, 1, figsize=(10, 6), sharex=True)
    figure_name = '.'.join([tr1.stats.network, tr1.stats.station, tr1.stats.channel[-1] + 
                 tr2.stats.channel[-1] + tr3.stats.channel[-1], figure_type])
    figure_path = os.path.join(figure_dir, figure_name)

    # plot the seismic record
    ax[0].plot(t, data1, color='black', label=tr1.stats.channel[-1])
    ax[0].set_xlim(0, np.max(t))
    ax[0].set_ylabel('Velocity (nm/s)', fontsize=14)
    ax[0].axvline(x=p_loc, color='blue', linestyle='--')
    ax[0].text(p_loc, np.min(data1), 'P', ha='right', color='blue', fontsize=15)
    ax[0].legend()  

    ax[1].plot(t, data2, color='black', label=tr2.stats.channel[-1])
    ax[1].set_xlim(0, np.max(t))
    ax[1].set_ylabel('Velocity (nm/s)', fontsize=14)
    ax[1].axvline(x=p_loc, color='blue', linestyle='--')
    ax[1].text(p_loc, np.min(data2), 'P', ha='right', color='blue', fontsize=15)
    ax[1].legend()  

    ax[2].plot(t, data3, color='black', label=tr3.stats.channel[-1])
    ax[2].set_xlim(0, np.max(t))
    ax[2].set_xlabel('Time (s)', fontsize=14)
    ax[2].set_ylabel('Velocity (nm/s)', fontsize=14)
    ax[2].axvline(x=p_loc, color='blue', linestyle='--')
    ax[2].text(p_loc, np.min(data3), 'P', ha='right', color='blue', fontsize=15)
    ax[2].legend()  

    try:
        s_loc = tr1.stats.sac.t2 - tr1.stats.sac.b
        ax[0].axvline(x=s_loc, color='red', linestyle='--')
        ax[0].text(s_loc, np.min(data1), 'S', ha='right', color='red', fontsize=15)
        ax[1].axvline(x=s_loc, color='red', linestyle='--')
        ax[1].text(s_loc, np.min(data2), 'S', ha='right', color='red', fontsize=15) 
        ax[2].axvline(x=s_loc, color='red', linestyle='--')
        ax[2].text(s_loc, np.min(data3), 'S', ha='right', color='red', fontsize=15)     

    except:
        pass

    plt.tight_layout()
    plt.savefig(figure_path)
    plt.close()

    return

def plot_waveform_list(waveform_list, figure_path, title, fig_size, ev_o, filter_range):
    """
    Plot waveforms from a list of waveform files with annotations for P and S arrivals.
    Args:
        - waveform_list (list): List of waveform file paths.
        - figure_path (str): Path to save the plotted figure.
        - title (str): Title of the plot.
        - fig_size (tuple): Figure size (width, height).
        - ev_o (UTCDateTime): Origin time of the event.
        - filter_range (list): Range for bandpass filter.
    Returns:
        - None
    """
    st = Stream()
    for waveform in waveform_list:
        st += read(waveform)

    if filter_range != None:
        st.filter('bandpass', freqmin=filter_range[0], freqmax=filter_range[1])

    # Get the minimum start time
    min_start_time = st[0].stats.starttime
    for tr in st:
        tr.stats.distance = tr.stats.sac['dist'] * 1000  # Convert distance to meters
        if utc2int(tr.stats.starttime) <= utc2int(min_start_time):
            min_start_time = tr.stats.starttime

    # Plot the waveform
    fig = plt.figure(figsize=fig_size)
    plt.rcParams['xtick.labelsize'] = 16  # Set the font size for tick labels on the x-axis
    plt.rcParams['ytick.labelsize'] = 16
    st.plot(type='section', orientation='horizontal', linewidth=1.0, show=False, fig=fig)

    # Label the station, P arrival, and S arrival
    ax = fig.axes[0]
    ax.set_ylabel('Distance (km)', fontsize=18)
    ax.set_xlabel('Time (s)', fontsize=18)

    for i in range(len(st)):

        tr = st[i]
        # Label the P arrivals
        p_pick = ev_o + tr.stats.sac['t1'] - min_start_time
        ax.text(p_pick, tr.stats.distance / 1e3, '|', ha='center', color='blue', fontsize=20)
        # Label the station name
        ax.text(0, tr.stats.distance / 1e3, tr.stats.sac['kstnm'], color='red', fontsize=14)

        # Label the S arrivals if exist
        try:
            s_pick = ev_o + tr.stats.sac['t2'] - min_start_time
            ax.text(s_pick, tr.stats.distance / 1e3, '|', ha='center', color='red', fontsize=20)
        except:
            pass

    plt.title(title, fontsize=20)
    plt.savefig(figure_path, dpi=300)
    plt.close()

    return

def plot_stft(waveform_path, filter_range, figure_dir, figure_type):
    """
    Plot seismograms in time domain and time-frequency domain.
    Args:
        - waveform_path (str): Path to the waveform file.
        - filter_range (list): Range for bandpass filter.
        - figure_dir (str): Directory to save the figure.
        - figure_type (str): Type of the figure file (e.g., 'png', 'jpg').
    Returns:
        - None
    """
    plt.rcParams['xtick.labelsize'] = 13  # Set the font size for tick labels on the x-axis
    plt.rcParams['ytick.labelsize'] = 13

    # read the data
    tr = read(waveform_path)[0]

    if filter_range != None:
        tr.filter('bandpass', freqmin=filter_range[0], freqmax=filter_range[1])

    # Extract data and sampling rate
    data = tr.data
    sampling_rate = tr.stats.sampling_rate
    t = np.arange(0, len(data)) / sampling_rate

    # P arrival location on the waveform
    p_loc = tr.stats.sac.t1 - tr.stats.sac.b

    # Plot the figure
    fig, ax = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    figure_name = '.'.join([tr.stats.network, tr.stats.station, tr.stats.channel, 'stft', figure_type])
    figure_path = os.path.join(figure_dir, figure_name)

    # plot the seismic record
    ax[0].plot(t, data, color='black')
    ax[0].set_xlim(0, np.max(t))
    ax[0].set_ylabel('Velocity (nm/s)', fontsize=15)
    ax[0].axvline(x=p_loc, color='blue', linestyle='--')
    ax[0].text(p_loc, np.min(data), 'P', ha='right', color='blue', fontsize=15)

    # Label the S arrivals if exist
    try:
        s_loc = tr.stats.sac.t2 - tr.stats.sac.b
        ax[0].axvline(x=s_loc, color='red', linestyle='--')
        ax[0].text(s_loc, np.min(data), 'S', ha='right', color='red', fontsize=15)
    except:
        pass

    # plot the spectrogram
    spectrogram(data=data, samp_rate=sampling_rate, axes=ax[1], cmap='hot')
    ax[1].set_ylabel('Frequency (Hz)', fontsize=15)
    ax[1].set_xlabel('Time (s)', fontsize=15)

    # Save the figure
    plt.tight_layout()
    plt.savefig(figure_path)
    plt.close()

    return

if __name__ == '__main__':

    plot_waveform_list(waveform_list=glob.glob('./event_data/*'), figure_path="./figures/waveform_plot.png", 
                     title="Waveform Plot", fig_size= (12, 8), ev_o= UTCDateTime(2021, 9, 17, 10, 10, 23, 13), filter_range=[0.1, 1])

    plot_3chns_event(filepath_chn1='./ne_rt_data/AEJ.1001.EHE', filepath_chn2='./ne_rt_data/AEJ.1001.EHN', 
                filepath_chn3='./ne_rt_data/AEJ.1001.EHZ', filter_range=[0.1, 1], figure_dir='./figures', figure_type='png')

    plot_3chns_event(filepath_chn1='./ne_rt_data/AEJ.1001.EHR', filepath_chn2='./ne_rt_data/AEJ.1001.EHT', 
                filepath_chn3='./ne_rt_data/AEJ.1001.EHZ', filter_range=[0.1, 1], figure_dir='./figures', figure_type='png')

    plot_stft(waveform_path='./ne_rt_data/AEJ.1001.EHZ', filter_range=[0.1, 1], figure_dir='./figures', figure_type='png')


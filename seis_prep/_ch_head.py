from obspy import UTCDateTime
from obspy.core import Trace
from obspy.io.sac.util import obspy_to_sac_header

def ch_head_event(trace, sta_info, ev_info, arrs):
    """
    Modify the header of an ObsPy trace object with event and station information.
    Args:
        - trace (Trace): The input ObsPy Trace object.
        - sta_info (list): List containing station information [net, sta, stlo, stla, stel].
        - ev_info (list): List containing event information [ev_o, evlo, evla, evdp, mag].
        - arrs (list): List of arrival times [p_arr, s_arr (optional)].
    Returns:
        - Trace: The modified ObsPy Trace object.
    """

    tr = trace.copy()
    # Add SAC header to the trace
    tr.stats.sac = obspy_to_sac_header(tr.stats)
    # print(arrs)
    # Extract and validate station information
    try:
        net, sta, stlo, stla, stel = sta_info[0], sta_info[1], \
                                     sta_info[2], sta_info[3], sta_info[4]
    except (IndexError, ValueError) as e:
        raise ValueError("Invalid station information provided.") from e

    # Extract and validate event information
    try:
        ev_o, evlo, evla, evdp, mag = ev_info[0], ev_info[1], \
                                      ev_info[2], ev_info[3], ev_info[4]
    except (IndexError, ValueError) as e:
        raise ValueError("Invalid event information provided.") from e

    # Extract and validate arrival times
    try:
        if len(arrs) == 2:
            p_arr, s_arr = arrs[0], arrs[1]
        elif len(arrs) == 1:
            p_arr = arrs[0]
            s_arr = None
        else:
            raise ValueError("Invalid number of arrival times provided.")
    except IndexError as e:
        raise ValueError("Invalid arrival times provided.") from e

    # Modify reference time in the SAC header
    event_time = UTCDateTime(ev_o)
    tr.stats.sac['nzyear'] = event_time.year
    tr.stats.sac['nzjday'] = event_time.julday
    tr.stats.sac['nzhour'] = event_time.hour
    tr.stats.sac['nzmin']  = event_time.minute
    tr.stats.sac['nzsec']  = event_time.second
    tr.stats.sac['nzmsec'] = int(event_time.microsecond / 1000)

    # Update station information in the SAC header
    tr.stats.sac.knetwk = net
    tr.stats.sac.kstnm = sta
    tr.stats.sac.stla = stla
    tr.stats.sac.stlo = stlo
    tr.stats.sac.stel = stel

    # Update event information in the SAC header
    tr.stats.sac.evla = evla
    tr.stats.sac.evlo = evlo
    tr.stats.sac.evdp = evdp
    tr.stats.sac.mag = mag

    # Update arrival times in the SAC header
    tr.stats.sac.t1 = p_arr
    if s_arr is not None:
        tr.stats.sac.t2 = s_arr

    tr.stats.sac.lcalda = True

    return tr

import os
from obspy import read
from obspy.signal.rotate import rotate_ne_rt

def rotate_ne2rt(filepath_n, filepath_e, filepath_r, filepath_t):
    """
    Rotate N and E channel to R and T channel.
    Args:
        - filepath_n (str): Path to the file containing the North component data.
        - filepath_e (str): Path to the file containing the East component data.
        - filepath_r (str): Path where the Radial component data will be saved.
        - filepath_t (str): Path where the Transverse component data will be saved.
    Returns:
        - None
    """
    try:
        # Read the N and E traces
        tr_n = read(filepath_n)[0]
        tr_e = read(filepath_e)[0]

        chn_prefix = tr_n.stats.channel[:2]

        # Copy the N and E traces to create R and T traces
        tr_r = tr_n.copy()
        tr_t = tr_e.copy()

        # Read the backazimuth and data
        baz = tr_n.stats.sac.baz
        data_n, data_e = tr_n.data, tr_e.data

        # Rotate the data
        data_r, data_t = rotate_ne_rt(n=data_n, e=data_e, ba=baz)

        # Assign the rotated data to radial and tangential traces
        tr_r.data = data_r
        tr_r.stats.channel = chn_prefix + 'R'
        tr_r.write(filepath_r, format='SAC')

        tr_t.data = data_t
        tr_t.stats.channel = chn_prefix + 'T'
        tr_t.write(filepath_t, format='SAC')

        print(f"Rotation successful: {filepath_r}, {filepath_t}")
        
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == '__main__':
    
    filepath_n = './ne_rt_data/AEJ.1001.EHN'
    filepath_e = './ne_rt_data/AEJ.1001.EHE'
    filepath_r = './ne_rt_data/AEJ.1001.EHR'
    filepath_t = './ne_rt_data/AEJ.1001.EHT'

    rotate_ne2rt(filepath_n=filepath_n, filepath_e=filepath_e, filepath_r=filepath_r, filepath_t=filepath_t)

import numpy as np
from obspy import read
from obspy.geodetics import gps2dist_azimuth

def rotate_ne2rt(filepath_nn, filepath_ee, filepath_ne, filepath_en, filepath_rr, filepath_tt):
    """
    Rotate NN, NE, EN, EE cross-correlations to get RR and TT cross-correlations.
    Refer to Fan-chi Lin et al., 2008, formula (1).
    Args:
        - filepath_nn (str): Path to the file of NN cross-correlation.
        - filepath_ee (str): Path to the file of EE cross-correlation.
        - filepath_ne (str): Path to the file of NE cross-correlation.
        - filepath_en (str): Path to the file of EN cross-correlation.
        - filepath_rr (str): Path to save the RR cross-correlation.
        - filepath_tt (str): Path to save the TT cross-correlation.
    Returns:
        - None
    """
    # Read the data
    tr_nn = read(filepath_nn)[0]
    tr_ee = read(filepath_ee)[0]
    tr_ne = read(filepath_ne)[0]
    tr_en = read(filepath_en)[0]

    # Read the location of stations, and calculate the azimuth and back azimuth
    evlo, evla, stlo, stla = tr_nn.stats.sac.evlo, tr_nn.stats.sac.evla, tr_nn.stats.sac.stlo, tr_nn.stats.sac.stla
    _, az, baz = gps2dist_azimuth(lat1=evla, lon1=evlo, lat2=stla, lon2=stlo, a=6378137.0, f=0.0033528106647474805)

    # Convert azimuths from degrees to radians
    theta = np.deg2rad(az)
    psi = np.deg2rad(baz)

    # Define rotation matrix as per the formula
    R = np.array([
        [-np.cos(theta) * np.cos(psi),  np.cos(theta) * np.sin(psi), -np.sin(theta) * np.sin(psi),  np.sin(theta) * np.cos(psi)],
        [-np.sin(theta) * np.sin(psi), -np.sin(theta) * np.cos(psi), -np.cos(theta) * np.cos(psi), -np.cos(theta) * np.sin(psi)],
        [-np.cos(theta) * np.sin(psi), -np.cos(theta) * np.cos(psi),  np.sin(theta) * np.cos(psi),  np.sin(theta) * np.sin(psi)],
        [-np.sin(theta) * np.cos(psi),  np.sin(theta) * np.sin(psi),  np.cos(theta) * np.sin(psi), -np.cos(theta) * np.cos(psi)]
    ])

    # Combine the data into a vector
    data_vector = np.array([
        tr_ee.data,
        tr_en.data,
        tr_nn.data,
        tr_ne.data
    ])

    # Perform the matrix multiplication to get the rotated data
    rotated_data = R @ data_vector

    # Assign the rotated data to RR and TT traces
    tr_rr = tr_nn.copy()
    tr_tt = tr_nn.copy()

    tr_tt.data = rotated_data[0]  # TT component
    tr_rr.data = rotated_data[1]  # RR component

    # Save the rotated traces
    tr_rr.write(filepath_rr, format='SAC')
    tr_tt.write(filepath_tt, format='SAC')

    return

if __name__ == '__main__':

    filepath_nn = './data_cc/X1.4501_X1.4502_NN_ls.SAC'
    filepath_ee = './data_cc/X1.4501_X1.4502_EE_ls.SAC'
    filepath_ne = './data_cc/X1.4501_X1.4502_NE_ls.SAC'
    filepath_en = './data_cc/X1.4501_X1.4502_EN_ls.SAC'

    filepath_rr = './data_cc/X1.4501_X1.4502_RR_ls.SAC'
    filepath_tt = './data_cc/X1.4501_X1.4502_TT_ls.SAC'

    rotate_ne2rt(filepath_nn=filepath_nn, filepath_ee=filepath_ee, filepath_ne=filepath_ne, filepath_en=filepath_en,
                 filepath_rr=filepath_rr, filepath_tt=filepath_tt)

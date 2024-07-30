from obspy.taup import TauPyModel

def cal_p_first_arrival(vel_model, dist_deg, evdp):
    """
    Calculate the first arrival time of P waves given a velocity model, distance, and event depth.
    Args:
        - vel_model (str): Name of the velocity model to use (e.g., 'iasp91', 'ak135').
        - dist_deg (float): Distance in degrees between the event and the station.
        - evdp (float): Depth of the event in kilometers.
    Returns:
        - float: The first arrival time of P waves in seconds, rounded to two decimal places.
    """
    model = TauPyModel(model=vel_model)

    # Get the P wave arrivals
    p_arrivals = model.get_travel_times(source_depth_in_km=evdp, 
                                        distance_in_degree=dist_deg, 
                                        phase_list=['P', 'Pn', 'Pg', 'p', 'Pdiff'])

    # Extract the arrival times
    p_arr_list = [arrival.time for arrival in p_arrivals]

    # Check if there are any arrivals
    if not p_arr_list:
        raise ValueError("No P wave arrivals found for the given parameters.")

    # Get the earliest arrival time
    p_arr = min(p_arr_list)

    # Return the earliest arrival time rounded to two decimal places
    return round(p_arr, 2)

def cal_s_first_arrival(vel_model, dist_deg, evdp):
    """
    Calculate the first arrival time of S waves given a velocity model, distance, and event depth.
    Args:
        - vel_model (str): Name of the velocity model to use (e.g., 'iasp91', 'ak135').
        - dist_deg (float): Distance in degrees between the event and the station.
        - evdp (float): Depth of the event in kilometers.
    Returns:
        - float: The first arrival time of S waves in seconds, rounded to two decimal places.
    """
    model = TauPyModel(model=vel_model)

    # Get the S wave arrivals
    s_arrivals = model.get_travel_times(source_depth_in_km=evdp, 
                                        distance_in_degree=dist_deg, 
                                        phase_list=['S', 'Sn', 'Sg', 's', 'Sdiff'])

    # Extract the arrival times
    s_arr_list = [arrival.time for arrival in s_arrivals]

    # Check if there are any arrivals
    if not s_arr_list:
        raise ValueError("No S wave arrivals found for the given parameters.")

    # Get the earliest arrival time
    s_arr = min(s_arr_list)

    # Return the earliest arrival time rounded to two decimal places
    return round(s_arr, 2)


if __name__ == '__main__':

    vel_model = 'iasp91'
    dist_deg = 30.0
    evdp = 10.0

    try:
        p_first_arrival = cal_p_first_arrival(vel_model, dist_deg, evdp)
        print(f"First P wave arrival time: {p_first_arrival} seconds")
    except ValueError as e:
        print(e)

    try:
        s_first_arrival = cal_s_first_arrival(vel_model, dist_deg, evdp)
        print(f"First S wave arrival time: {s_first_arrival} seconds")
    except ValueError as e:
        print(e)

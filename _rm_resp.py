import os
from obspy import UTCDateTime
from obspy import read
from obspy import read_inventory
from obspy.io.xseed.core import _read_resp
from obspy.io.xseed import Parser
from obspy.io.sac.sacpz import _write_sacpz
import subprocess

def xml2pz(xml_path, pz_dir):
    """
    Converts an XML inventory file to individual PZ (Pole-Zero) files for each channel.
    Args:
        - xml_path (str): Path to the XML file containing the inventory.
        - pz_dir (str): Directory where the PZ files will be saved.
    Returns:
        - None
    """
    # Read the inventory from the XML file
    inventory = read_inventory(xml_path)

    # Iterate over networks, stations, and channels
    for network in inventory:
        for station in network:
            for channel in station:

                # print(network, station, channel)
                new_inventory = inventory.select(network=network.code, station=station.code, location='00', channel=channel.code)
                if len(new_inventory) == 0:
                    new_inventory = inventory.select(network=network.code, station=station.code, channel=channel.code)

                # Define the output PZ file name
                pz_filename = f"{network.code}.{station.code}.{channel.code}.PZ"
                pz_file_path = os.path.join(pz_dir, pz_filename)
                
                _write_sacpz(new_inventory, pz_file_path)

# # some problems need to solve
# def xml2resp(xml_path, resp_dir):
#     """
#     Converts an XML inventory file to individual resp files for each channel.
#     Args:
#         - xml_path (str): Path to the XML file containing the inventory.
#         - resp_dir (str): Directory where the PZ files will be saved.
#     Returns:
#         - None
#     """

#     inventory = read_inventory(xml_path)
#     datetime = UTCDateTime("2005-08-24T00:00:00")

#     for network in inventory:
#         for station in network:
#             for channel in station:

#                 chn_name = f"{network.code}.{station.code}.{channel.location_code}.{channel.code}"
#                 print(chn_name)
#                 response = inventory.get_response(chn_name, datetime)
#                 resp_filename = f"{resp_dir}/RESP.{chn_name}"

#                 with open(resp_filename, 'w') as resp_file:
#                     resp_file.write(response)

def resp2pz(resp_path, pz_dir):
    """
    Converts an resp inventory file to individual PZ (Pole-Zero) files for each channel.
    Args:
        - resp_path (str): Path to the resp file containing the inventory.
        - pz_dir (str): Directory where the PZ files will be saved.
    Returns:
        - None
    """
    inventory = _read_resp(resp_path)

    # Iterate over networks, stations, and channels
    for network in inventory:
        for station in network:
            for channel in station:

                new_inventory = inventory.select(network=network.code, station=station.code, location='00', channel=channel.code)
                if len(new_inventory) == 0:
                    new_inventory = inventory.select(network=network.code, station=station.code, channel=channel.code)
                    
                pz_filename = f"{network.code}.{station.code}.{channel.code}.PZ"
                pz_file_path = os.path.join(pz_dir, pz_filename)

                _write_sacpz(new_inventory, pz_file_path)         


def rm_response_resp(file_path, out_file_path, resp_path, filter_range):
    """
    Removes the instrument response from a seismic data file using an RESP file.
    Args:
        - file_path (str): Path to the seismic data file.
        - out_file_path (str): Path for output file.
        - resp_path (str): Path to the RESP file.
        - filter_range (str): Frequency range to apply (e.g., '0.01 0.02 10 20').
    """
    # Set environment variable to suppress SAC copyright display
    os.putenv("SAC_DISPLAY_COPYRIGHT", '0')
    # Create the SAC command script
    sac_commands = f"""
    r {file_path}
    rmean
    rtrend
    taper
    trans from evalresp fname {resp_path} to vel freq {filter_range}
    w {out_file_path}
    q
    """
    # Run the SAC command script
    subprocess.Popen(['sac'], stdin=subprocess.PIPE).communicate(sac_commands.encode())

def rm_response_pz(file_path, out_file_path, pz_path, filter_range):
    """
    Removes the instrument response from a seismic data file using a pole-zero file.
    Args:
        - file_path (str): Path to the seismic data file.
        - out_file_path (str): Path for output file.
        - pz_path (str): Path to the pole-zero file.
        - filter_range (str): Frequency range to apply (e.g., '0.01 0.02 10 20').
    """
    # Set environment variable to suppress SAC copyright display
    os.putenv("SAC_DISPLAY_COPYRIGHT", '0')
    # Create the SAC command script
    sac_commands = f"""
    r {file_path}
    rmean
    rtrend
    taper
    trans from polezero subtype {pz_path} to vel freq {filter_range}
    mul 1.0e9
    w {out_file_path}
    q
    """
    # Run the SAC command script
    subprocess.Popen(['sac'], stdin=subprocess.PIPE).communicate(sac_commands.encode())

if __name__ == '__main__':

    import matplotlib.pyplot as plt

    ####################### 1 RESP VS PZ for chn stations #######################

    # sac_file1 = './resp_pz/SD.TIA.BHZ'
    # resp_file = './resp_pz/RESP.SD.TIA.00.BHZ'
    # pz_file = './resp_pz/SD.TIA.BHZ.PZ'

    # resp2pz(resp_path=resp_file, pz_dir='./resp_pz')

    # rm_response_resp(file_path=sac_file1, out_file_path='./resp_pz/resp_SD.TIA.BHZ', resp_path=resp_file, filter_range='0.005 0.01 40 45')
    # rm_response_pz(file_path=sac_file1, out_file_path='./resp_pz/pz_SD.TIA.BHZ', pz_path=pz_file, filter_range='0.005 0.01 40 45')

    # tr_raw = read('./resp_pz/SD.TIA.BHZ')[0]
    # tr_pz = read('./resp_pz/pz_SD.TIA.BHZ')[0]
    # tr_resp = read('./resp_pz/resp_SD.TIA.BHZ')[0]

    # plt.plot(tr_raw.data, color='black', label='raw')
    # plt.plot(tr_pz.data, color='red', label='pz')
    # plt.plot(tr_resp.data, color='blue', label='resp')

    # plt.legend()
    # plt.show()


    # ####################### 2 RESP VS PZ for iris stations #######################

    sac_file3 = './resp_pz_iris/IU.COLA.BHZ.SAC'
    resp_file = './resp_pz_iris/RESP.IU.COLA.00.BHZ'
    pz_file = './resp_pz_iris/IU.COLA.BHZ.PZ'

    resp2pz(resp_path=resp_file, pz_dir='./resp_pz_iris')

    rm_response_resp(file_path=sac_file3, out_file_path='./resp_pz_iris/resp_IU.COLA.BHZ.SAC', resp_path=resp_file, filter_range='0.01 0.05 5 7')
    rm_response_pz(file_path=sac_file3, out_file_path='./resp_pz_iris/pz_IU.COLA.BHZ.SAC', pz_path=pz_file, filter_range='0.01 0.05 5 7')

    tr_raw = read('./resp_pz_iris/IU.COLA.BHZ.SAC')[0]
    tr_raw.detrend('demean')

    tr_pz = read('./resp_pz_iris/pz_IU.COLA.BHZ.SAC')[0]
    tr_resp = read('./resp_pz_iris/resp_IU.COLA.BHZ.SAC')[0]

    plt.plot(tr_raw.data, color='black', label='raw')
    plt.plot(tr_pz.data, color='red', label='pz')
    plt.plot(tr_resp.data, color='blue', label='resp')

    # print(tr_pz.data == tr_resp.data)

    plt.legend()
    plt.show()


    ####################### 3 XML VS PZ for iris stations #######################
    # sac_file2 = './xml_pz/IC.LSA.BHZ.SAC'
    # xml_file = './xml_pz/IC.LSA.xml'
    # pz_file2 = './xml_pz/IC.LSA.BHZ.PZ'

    # xml2pz(xml_path=xml_file, pz_dir='./xml_pz')
    # rm_response_pz(file_path=sac_file2, out_file_path='./xml_pz/pz_IC.LSA.BHZ.SAC', pz_path=pz_file2, filter_range='0.005 0.01 40 45')

    # st = read(sac_file2)
    # tr_raw = st[0]

    # inventory = read_inventory(xml_file)
    # st_xml = st.remove_response(inventory=inventory, water_level=None, pre_filt = [0.005, 0.01, 40, 45])
    # st_xml.write('./xml_pz/xml_IC.LSA.BHZ.SAC', format='SAC')

    # tr_xml = read('./xml_pz/xml_IC.LSA.BHZ.SAC')[0]
    # tr_pz = read('./xml_pz/pz_IC.LSA.BHZ.SAC')[0]

    # plt.plot(tr_raw.data, color='black', label='raw')
    # plt.plot(tr_pz.data / 1e9, color='red', label='pz')
    # plt.plot(tr_xml.data, color='blue', label='xml')
    # # plt.plot(tr_xml.data - tr_pz.data/1e9, color='black', label='xml_pz')
    # # plt.plot(tr_pz.data/1e9 - tr_raw.data, color='red', label='pz_raw')
    # # plt.plot(tr_xml.data - tr_raw.data, color='blue', label='xml_raw')


    # print(tr_raw.data == tr_xml.data)

    # plt.legend()
    # plt.show()    

    ####################### 4 XML VS resp for iris stations #######################
    # sac_file4 = './xml_pz/IC.LSA.BHZ.SAC'
    # xml_file = './xml_pz/IC.LSA.xml'
    # resp_file = './xml_pz/IC.LSA.BHZ.PZ'

    # xml2resp(xml_path=xml_file, resp_dir='./xml_pz')
    # rm_response_resp(file_path=sac_file4, out_file_path='./xml_pz/resp_IC.LSA.BHZ.SAC', resp_path=resp_file, filter_range='0.005 0.01 40 45')

    # st = read(sac_file4)
    # tr_raw = st[0]

    # tr_xml = st.remove_response(read_inventory(xml_file), water_level=None, pre_filt = [0.005, 0.01, 40, 45])[0]
    # tr_resp = read('./xml_pz/resp_IC.LSA.BHZ.SAC')[0]

    # plt.plot(tr_resp.data / 1e9, color='red', label='pz')
    # plt.plot(tr_xml.data, color='blue', label='xml')
    # plt.plot(tr_raw.data, color='black', label='raw')

    # print(tr_resp.data == tr_xml.data)

    # plt.legend()
    # plt.show()    


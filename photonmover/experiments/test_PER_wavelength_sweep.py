#%%
# Interfaces/instruments necessary for the experiment
# - You use an Interface if any instrument of that category can be used
# - You use a specific instrument if you can only use that specific model
from photonmover.instruments.Power_meters.Thorlabs import ThorlabsPowerMeter
# from photonmover.instruments.Power_Supplies.KeysightE36106B import KeysightE36106B
from photonmover.instruments.Source_meters.Keithley2635A import Keithley2635A
from photonmover.instruments.Power_Supplies.AgilentE3633A import AgilentE3633A

# General imports
import time
import numpy as np
import csv
from scipy import io
import os
import matplotlib.pyplot as plt

from ctypes import create_string_buffer, c_int16, c_double, byref, c_char_p, c_bool
#%%
# Power meter sensors
XPOL_SENSOR = b'S122C'  # 800-1700nm
PPOL_SENSOR = b'S122C'  # 800-1700nm

def load_wavvolt(filename):
    """
    Read wavelength over voltage data from mat file
    """
    full_filename = filename + ".mat"
    a = io.loadmat(full_filename)
    wav_list = a['wav_select'][0]
    volt_list = a['volt_select'][0]
    return [volt_list, wav_list]

def load_delayvolt(delay_dir,filename):
    """
    Read wavelength over voltage data from mat file
    """
    full_filename = delay_dir + "\\" + filename + ".mat"
    a = io.loadmat(full_filename)
    VOA_interp = a['VOA_est_interp'][0]
    VCSEL_interp = a['V_MEMSinterp'][0]
    return [VOA_interp, VCSEL_interp]

def interp_VOA(VCSEL_set, VOA_interp, VCSEL_interp):
    VOA_set = []
    for i, volt in enumerate(VCSEL_set):
        ind = np.argmin(np.abs(VCSEL_interp - VCSEL_set[i]))
        VOA_set.append(VOA_interp[ind])
    return  np.array(VOA_set)


def perform_experiment(instrument_dict, params, filename=None):
    """
    Performs the experiment, and saves the relevant data (if there is any)
    to the specified file (if given)
    :param params: dictionary of the parameters necessary for the experiment.
    :param filename: if specified, the data is saved in the specified file.
    :return: [meas_volt_list, PER]
    """

    voa = instrument_dict["voa"]
    ps = instrument_dict["ps"]
    pm1 = instrument_dict["pm1"]
    pm2 = instrument_dict["pm2"]
    
    [pm1_p0, _] = pm1.get_powers()
    [pm2_p0, _] = pm2.get_powers()
    
    if pm1_p0 > pm2_p0:
        pm_ppol = pm1
        pm_xpol = pm2
    elif pm2_p0 < pm1_p0:
        pm_ppol = pm2
        pm_xpol = pm1
    else: #shouldn't happen, but randomly assign if does
        pm_ppol = pm1
        pm_xpol = pm2
    
    VCSEL_wavelength = params["VCSEL_wavelength"]
    voltage = params["voltage"]
    if voa is not None:
        voa_set = params["voa_set"]

    # print('Measured VCSEL pm wavelength is %d nm' % meas_xpol_wavelength.value)
    ppol_power_list = []
    xpol_power_list = []
    meas_volt_list = []

    # Sweep VOA voltage, set power meter wavelength and get power
    #Set first voltage and allow time to settle
    ps.set_voltage(voltage[0])
    if voa is not None:
        voa.set_voltage(voa_set[0])
    time.sleep(1)

    for (ind, volt) in enumerate(voltage):

        # Set power meter wavelength
        pm_ppol.setWavelength(c_double(VCSEL_wavelength[ind]))
        pm_xpol.setWavelength(c_double(VCSEL_wavelength[ind]))

        # Wait s
        time.sleep(0.5)

        print('Setting power supply voltage to %.4f V...' % volt)
        # Set the voltage
        ps.set_voltage(volt)
        # if params["voa_set"] is not None:
        if voa is not None:
            voa.set_voltage(voa_set[ind])
        time.sleep(1)

        meas_volt = ps.measure_voltage()
        print('Voltage set to %0.4f V' % meas_volt)
        meas_volt_list.append(meas_volt) #[V]

        # Wait s
        time.sleep(2.0)

        # Get VCSEL power at cross and parallel polarizations
        [ppol_power, _] = pm_ppol.get_powers()
        ppol_power_list.append(ppol_power) #[W]
        print('Measured PPol power is %.6f mW' % (ppol_power*1e3))

        [xpol_power, _] = pm_xpol.get_powers()
        xpol_power_list.append(xpol_power)  # [W]
        print('Measured XPol power is %.6f mW' % (xpol_power * 1e3))

    PER = 10*np.log10(np.array(ppol_power_list) / np.array(xpol_power_list))

    print('Finished PER wavelength sweep')
    print('-----------------------------')

    if filename is not None:
        # Save the data in a csv file
        time_tuple = time.localtime()
        complete_filename = "%s-%d-%d-%d.csv" % (filename,
                                                        time_tuple[0],
                                                        time_tuple[1],
                                                        time_tuple[2])

        with open(complete_filename, 'w+') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(meas_volt_list)  #[V]
            writer.writerow(VCSEL_wavelength) #[nm]
            writer.writerow(ppol_power_list)  #[W]
            writer.writerow(xpol_power_list)  #[W]
            writer.writerow(PER) #[dB]
            if voa is not None:
                writer.writerow(voa_set)  #[V]

    return [meas_volt_list, PER]



#%%
# ------------------------------------------------------------
# SAFETY LIMITS
i_limit = 10e-9 #current limit

# POWER METER SETTINGS
# Check that power meter sensors (declared as global variables up top) are accurate
# Note: pump sensor and VCSEL sensor can be the same type
pump_wavelength = 1076 #nm

#OTHER DEVICE PARAMETERS
wavvolt_dir = os.path.join(os.getcwd(), "wavvolt")
delayvolt_dir = os.path.join(os.getcwd(), "delayvolt")

device ='dev1b_SS_delayvolt16_CW159mA_1'

# EXPERIMENT PARAMETERS
vary_delay = 1
pump_power = 0#14 # [mW]
IL = 0.5#0.50

wavvolt_file = "wavvolt_dev1b_15C_OE1076_0.00mW_2025-11-1_delayvolt16"
delayvolt_file = "delayvolt16"
fwname = os.path.join(wavvolt_dir, wavvolt_file)


[volt_list, VCSEL_wavelength] = load_wavvolt(fwname)
if vary_delay:
    [VOA_interp, VCSEL_interp] = load_delayvolt(delayvolt_dir, delayvolt_file)
    voa_set = interp_VOA(volt_list, VOA_interp, VCSEL_interp)
    voa = AgilentE3633A()
    voa.initialize()
else:
    voa_set = None

# VCSEL_wavelength = VCSEL_wavelength[volt_list < end_voltage]
# volt_list = volt_list[volt_list < end_voltage]
print(volt_list)
# ------------------------------------------------------------
#%%
# INSTRUMENTS
# ps = KeysightE36106B(current_limit=i_limit)
ps = Keithley2635A(current_compliance=i_limit, voltage_compliance=-76)
pm1 = ThorlabsPowerMeter()
pm2 = ThorlabsPowerMeter()

# Set up power supply
ps.initialize()

# Set up power meters
pm1.initialize(deviceIndex=0)
pm2.initialize(deviceIndex=1)
pm1.setPowerAutoRange(c_int16(1))  #enable autorange
pm2.setPowerAutoRange(c_int16(1)) #enable autorange

file_name = "PER_%s_%d-%dV_GS%dnm_%2.3fmW" % (device, volt_list[0], volt_list[-1], pump_wavelength, pump_power*IL) # Filename where to save csv data

# SET UP THE EXPERIMENT
instr_dict = {
    "pm1": pm1,
    "pm2": pm2,
    "ps": ps,
    "voa": voa
}
params = {"xpol_sensor": XPOL_SENSOR, "ppol_sensor": PPOL_SENSOR, "VCSEL_wavelength": VCSEL_wavelength * 1e9,
                "voltage": volt_list, "voa_set": voa_set}

# RUN IT
[volt_list, PER_list] = perform_experiment(instr_dict, params, filename=file_name)

# CLOSE INSTRUMENTS
pm1.close()
pm2.close()
ps.close()
voa.close()
#%% Plot data

fig, ax = plt.subplots(1,1, figsize=(6,4), tight_layout=True)
ax.plot(volt_list, PER_list)
ax.set_xlabel("MEMS Voltage (V)")
ax.set_ylabel("PER (dB)")
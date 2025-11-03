#%%
from photonmover.instruments.Oscilloscopes.HP54750A import HP54750A
from photonmover.instruments.Source_meters.Keithley2635A import Keithley2635A
from photonmover.instruments.Power_Supplies.AgilentE3633A import AgilentE3633A
from photonmover.instruments.Optical_spectrum_analyzers.HP70951B import HP70951B

# General imports
import time
import os
import numpy as np
import csv
from scipy import io
from scipy.signal import find_peaks
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt

#%%

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
    to the specified file (if given) - assumes osa set to desired display parameters prior to running experiment
    :param params: dictionary of the parameters necessary for the experiment.
    :param filename: if specified, the data is saved in the specified file.
    :return: [meas_volt_list, pk_wl_list]
    """
    ps = instrument_dict["ps"]
    voa = instrument_dict["voa"]
    det = instrument_dict["det"]
    
    voltage_list = params["voltage_list"]
    voa_set = params["voa_set"]

    meas_volt_list = []
    t_peak = []

    if isinstance(det, HP54750A):
        # Read trace from sampling oscilloscope
        full_filename = params["detector"] + "_" + filename
        [t, _] = det.read_waveform([1])
        trace_len = np.array(t).shape[1]
        trace_array = np.array(np.zeros((len(voltage_list), int(trace_len))))
        x_axis = t[0]
        t_sleep = 10
    elif isinstance(det, HP70951B):
        [wavelength, _] = det.read_data()
        trace_len = np.array(wavelength).shape[0]
        trace_array = np.array(np.zeros((len(voltage_list), int(trace_len))))
        x_axis = wavelength
        t_sleep = 0

    # Sweep power supply voltage and get osa trace
    for ind, volt in enumerate(voltage_list):
        print('Setting power supply to %.4f V...' % volt)
        # Set the voltage
        ps.set_voltage(volt)
        voa.set_voltage(voa_set[ind])

        meas_volt = ps.measure_voltage()
        print('Power Supply voltage set to %0.4f V' % meas_volt)
        meas_volt_list.append(meas_volt) #[V]

        # Wait [s]
        time.sleep(t_sleep)

        #Read osc trace and store in wavelength_array
        if isinstance(det, HP54750A):
            # Read trace from sampling oscilloscope
            full_filename = params["detector"] + "_" + filename
            [_, waveform] = det.read_waveform([1])
            trace_array[ind, :] = np.reshape(np.array(waveform), (1, int(trace_len)))

        elif isinstance(det, HP70951B):
            full_filename = params["detector"] + "_" + filename
            [_, waveform] = det.read_data()
            trace_array[ind, :] = np.reshape(np.array(waveform), (1, int(trace_len)))

        #Add code to measure distance between first two peaks
        win_len = 15
        poly_order = 3
        smooth_waveform = savgol_filter(np.array(waveform), win_len, poly_order)
        pk_ind = np.argmax(smooth_waveform)

        t_peak.append(np.array(x_axis)[pk_ind])
        # pks_ind, _ = find_peaks(smooth_waveform[0,:])
        # print(pks_ind)
        # delay = np.diff(t[pks_ind])
        # delay_list.append(delay)


    print(np.shape(trace_array))
    print('Finished voltage sweep')
    print('-----------------------------')


    if filename is not None:
        # Save the data in a csv file
        time_tuple = time.localtime()
        complete_filename = "%s-%d-%d-%d_%d-%d-%d.csv" % (filename,
                                                        time_tuple[0],
                                                        time_tuple[1],
                                                        time_tuple[2],
                                                        time_tuple[3],
                                                        time_tuple[4],
                                                        time_tuple[5])

        with open(complete_filename, 'w+') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(np.array(x_axis[:]))  #[s]
            for ind in range(len(voltage_list)):
                writer.writerow(trace_array[ind, :])  #[V]

                # Save the parameters in a .mat file
            time_tuple = time.localtime()
            params_filename = "params_%s_%s_%d-%d-%d.mat" % (
                params["detector"],
                filename,
                time_tuple[0],
                time_tuple[1],
                time_tuple[2])

            io.savemat(params_filename, {'voltage_list': voltage_list
                                            })
    return [meas_volt_list, t_peak]


        
#%%
# ------------------------------------------------------------
# SAFETY LIMITS
i_limit = 10e-9  # current limit

# OTHER PARAMETERS
detector = 'osc'
device = 'dev1b_delayvolt16_noBOA_04'
pump_laser = 'OE1076CW159mA' #'OEland1076' #'OEland1038' #'CW976'
pump_power = 0 #mW
IL = 0.62
RBW = 0.1 #nm
temp = 15 #C

# EXPERIMENT PARAMETERS
init_voltage = 0  # [V]
end_voltage = 60 # [V]
increment = 1  # Voltage increment
voltage_list = np.arange(init_voltage, end_voltage+increment, increment) #end_voltage+1 or will stop at end_voltage-1
voltage_list = np.append(voltage_list, [60.5, 61, 61.5, 61.9])
print(voltage_list)
delayvolt_file = "delayvolt16"
delay_dir = os.path.join(os.getcwd(),"delayvolt")
# ------------------------------------------------------------
#%%
# INSTRUMENTS
ps = Keithley2635A(current_compliance=10e-9, voltage_compliance=73) #A, V
if detector == "osc":
    det = HP54750A()
elif detector == "osa":
    det = HP70951B()
# osc = HP54750A()
voa = AgilentE3633A(current_limit=i_limit)

[VOA_interp, VCSEL_interp] = load_delayvolt(delay_dir, delayvolt_file)
voa_set = interp_VOA(voltage_list, VOA_interp, VCSEL_interp)

# Initialize instruments
ps.initialize()
det.initialize()
# osc.initialize()
voa.initialize()

# file_name = 'LL_dev1_50V_pm1258_Solstis980nm'  # Filename where to save csv data
file_name = "tuningDelay_%s_%s_%d-%dV_%dC_%s_%3.2fmW_RBW%3.2fnm" % (
    detector, device, init_voltage, end_voltage, temp, pump_laser, pump_power*IL, RBW)  # Filename where to save csv data

# SET UP THE EXPERIMENT
instr_dict = {
    "det": det,
    "ps": ps,
    "voa": voa
    }
 
params = {"detector": detector, "voltage_list": voltage_list, "voa_set": voa_set, "device": device, "pump_laser": pump_laser, "pump_power": pump_power, "temp": temp}

# RUN IT
[meas_volt_list, t_peak] = perform_experiment(instr_dict, params, filename=file_name)

# CLOSE INSTRUMENTS
det.close()
ps.close()

#%% Plot data
fig, ax = plt.subplots(1,1, figsize=(6,4), tight_layout=True)
ax.plot(meas_volt_list, t_peak)
ax.set_xlabel("MEMS Voltage (V)")
ax.set_ylabel("Peak time (ps)")

# %%

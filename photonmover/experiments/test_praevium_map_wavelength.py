#%% Lazy rewrite of praevium_map_wavelength.py for Yokogawa OSA using Gemini Pro
# Interfaces/instruments necessary for the experiment
# - You use an Interface if any instrument of that category can be used
# - You use a specific instrument if you can only use that specific model
# from photonmover.instruments.Optical_spectrum_analyzers.HP70951B import HP70951B
from photonmover.instruments.Source_meters.Keithley2635A import Keithley2635A

# General imports
import pyvisa
import time
import numpy as np
import csv
from scipy import io
import matplotlib.pyplot as plt


class YokogawaAQ6375:
    def __init__(self):
        self.osa = None
    
    def initialize(self, GPIB_ADDRESS='GPIB0::1::INSTR', timeout=50000):
        """
        Connects to Yokogawa AQ6375 given resource_name (e.g., 'GPIB0::1::INSTR')
        """
        try:
            rm = pyvisa.ResourceManager()
            self.osa = rm.open_resource(GPIB_ADDRESS, timeout=timeout)
            print('Opening connnection to Yokogawa AQ6375 OSA')
        except BaseException:
            raise ValueError('Cannot connect to the Yokogawa AQ6375 OSA')
        # Verify connection
        print(f"Connected to: {self.osa.query('*IDN?').strip()}")
       
        
    def close(self):
        self.osa.close()
            
    def read_data(self):
        """
        Performs a single sweep and returns wavelength and power data.
        tuple: (wavelengths, power_levels) as numpy arrays
        """
        try:
            # 1. Set to Single Sweep Mode and Trigger
            self.osa.write(":FORM:DATA ASCII")
            self.osa.write(":SENS:SWE:MODE SING") 
            self.osa.write("*CLS")                 # Clear status registers
            self.osa.write(":INIT:IMM")           # Trigger immediate sweep
            
            # 2. Wait for sweep completion
            # *OPC? returns '1' only after the sweep is finished
            print("Sweeping...")
            while True:
                try:
                    if self.osa.query("*OPC?").strip() == "1":
                        break
                except pyvisa.errors.VisaIOError:
                    # For very slow sweeps, the query might timeout; we keep waiting
                    time.sleep(1)
            
            # 3. Retrieve Data
            # :TRACe:X? TRA returns wavelength data in meters
            # :TRACe:Y? TRA returns power data (usually dBm)
            print("Transferring data...")
            raw_x = self.osa.query(":TRACe:X? TRA")
            raw_y = self.osa.query(":TRACe:Y? TRA")
            
            # print(f"raw_x: {raw_x}")
            raw_x = raw_x.replace("m","-") #Yokogawa randomly places m by some exponents
            # print(f"raw_y: {raw_y}")
            raw_y = raw_y.replace("/","-") # Y trace sometimes has '/' at begining
            
            # 4. Convert comma-separated strings to float arrays
            wavelengths = np.fromstring(raw_x, sep=',')
            powers = np.fromstring(raw_y, sep=',')
            
            return wavelengths, powers

        except Exception as e:
            print(f"Communication Error: {e}")
            return None, None
    
    def set_rbw(self, rbw):
        self.osa.write(f":SENSe:BANDwidth:RESolution {rbw}NM")
        
    def set_center_wavelength(self, wavelength_nm):
        self.osa.write(f":DISP:TRACe:X:CENTer {wavelength_nm}NM")
    
    def set_span(self, span):
        self.osa.write(f":SENSe:WAV:SPAN {span}NM")
        
    def set_ref_level(self, ref_level):
        self.osa.write(f":DISP:TRACe:Y1:RLEV {ref_level}DBM")
        
    def set_num_points(self, num_points):
        self.osa.write(f"SENSe:SWEep:POINts {num_points}")
        
    def get_rbw(self):
        rbw = self.osa.query_ascii_values(":SENS:BAND:RES?")
        return rbw[0]
    
    def get_center_wavelength(self):
        return self.osa.query_ascii_values(":DISP:TRACe:X:CENTer?")[0]
    
    def get_span(self):
        return self.osa.query_ascii_values(":SENSe:WAV:SPAN?")[0]
    
    def get_ref_level(self):
        return self.osa.query_ascii_values(":DISP:TRACe:Y1:RLEV?")[0]
    
    def get_num_points(self):
        return self.osa.query_ascii_values("SENSe:SWEep:POINts?")[0]
    
    def get_osa_parameters(self):
        """
        Read center wavelength, span, reference level, and RBW
        :return: [trace_len, osa_params] = dict(center_wavelength, span, reference_level, RBW, sensitivity)
        """
        center_wavelength = self.get_center_wavelength()
        span = self.get_span() #[m]
        reference_level = self.get_ref_level() #[dBm]
        RBW = self.get_rbw() #[m]
        # sensitivity = self.osa.query_ascii_values('SENS?') #[dBm]

        # osa_params = [center_wavelength, span, reference_level, RBW]
        osa_params = {
            'center wavelength':    center_wavelength,
            'span':                 span,
            'reference level':      reference_level,
            'rbw':                  RBW, 
            'sensitivity':          np.nan #For backward compatibility
        }

        trace_len = self.get_num_points()

        return [trace_len, osa_params]
        
def save_wavvolt(filename, volt, peak_wl, f_interp=5000):
        """
        Interpolates voltage vs peak wavelength curve, saves wavvolt file with variables:
        -volt_select: original measured voltage from experiment
        -wav_select: original measured wavelengths from experiment
        -volt_interp: interpolated voltage
        -wav_interp: interpolated wavelength
        """
        volt_interp = np.linspace(volt[0], volt[-1], f_interp)
        peak_interp = np.interp(volt_interp, volt, peak_wl)

        time_tuple = time.localtime()
        wavvolt_filename = "wavvolt_%s_%d-%d-%d.mat" % (
            filename,
            time_tuple[0],
            time_tuple[1],
            time_tuple[2])
        io.savemat(wavvolt_filename, {'volt_select': volt,
                                     'wav_select': peak_wl,
                                      'volt_interp': volt_interp,
                                      'peak_interp': peak_interp
                                     })
                            
def perform_experiment(params, instrument_dict: dict, filename=None):
        """
        Performs the experiment, and saves the relevant data (if there is any)
        to the specified file (if given) - assumes osa set to desired display parameters prior to running experiment
        :param params: dictionary of the parameters necessary for the experiment.
        :param filename: if specified, the data is saved in the specified file.
        :return: [meas_volt_list, pk_wl_list]
        """

        voltage_list = params["voltage_list"]
        osa = instrument_dict['osa']
        ps = instrument_dict['ps']

        meas_volt_list = []
        pk_wl_list = []
        [trace_len, osa_settings_list] = osa.get_osa_parameters()

        spectrum_array = np.array(np.zeros((len(voltage_list), int(trace_len))))
        [wavelength_list, _] = osa.read_data()

        # Sweep power supply voltage and get osa trace
        for ind, volt in enumerate(voltage_list):

            print('Setting power supply to %.4f V...' % volt)
            # Set the voltage
            ps.set_voltage(volt)

            meas_volt = ps.measure_voltage()
            print('Power Supply voltage set to %0.4f V' % meas_volt)
            meas_volt_list.append(meas_volt) #[V]

            # Wait [s]
            time.sleep(0.5)

            #Read osa spectrum and store in wavelength_array
            [_, amp] = osa.read_data()
            spectrum_array[ind, :] = np.reshape(np.array(amp), (1, int(trace_len)))
            pk_ind = np.argmax(spectrum_array[ind,:])
            pk_wl_list.append(wavelength_list[pk_ind])

        print(np.shape(spectrum_array))
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
                writer.writerow(wavelength_list)  #[m]
                for ind in range(len(voltage_list)):
                    writer.writerow(spectrum_array[ind, :])  #[W]

                    # Save the parameters in a .mat file
                time_tuple = time.localtime()
                params_filename = "params_%s_%d-%d-%d.mat" % (
                    filename,
                    time_tuple[0],
                    time_tuple[1],
                    time_tuple[2])

                io.savemat(params_filename, {'osa_settings_list': osa_settings_list,
                                             'voltage_list': voltage_list,
                                             'meas_volt_list': meas_volt_list
                                             })

                wavvolt_filename = "%s_%dC_%s_%3.2fmW" % (
                    params["device"],
                    params["temp"],
                    params["pump_laser"],
                    params["pump_power"]*params["IL"]
                )

                save_wavvolt(wavvolt_filename, np.array(meas_volt_list), np.array(pk_wl_list))

        return [meas_volt_list, pk_wl_list]
    
    
#%% Set parameters
# SAFETY LIMITS
i_limit = 10e-9  # current limit

# OTHER PARAMETERS
device = 'dev1b'
pump_laser = 'OE1076CW159mA_1.9V_BOA' #'OEland1038' #'CW976'
pump_power = 0 #10.2 #mW
IL = 1#0.62#0.52
temp = 15 #C

# OSA PARAMS
RBW = 0.1 #nm
start_wav = 1230 #nm
stop_wav = 1370 #nm
span = stop_wav - start_wav
center_wav = span/2 + start_wav
ref_level = -15.7 #dBm

# EXPERIMENT PARAMETERS
init_voltage = 0  # [V]
end_voltage = 69 # [V]
increment = 1  # Voltage increment
voltage_list = np.arange(init_voltage, end_voltage+increment, increment) #end_voltage+1 or will stop at end_voltage-1

#REMOVE IF END VOLTAGE CHANGES!
voltage_list = np.concatenate((voltage_list, [69.5, 70, 70.5, 71, 71.3, 71.5]))

print(f"voltage_list: {voltage_list}")
#%% Initialize instruements 
# ------------------------------------------------------------

# INSTRUMENTS
ps = Keithley2635A(current_compliance=i_limit, voltage_compliance=73) #A, V
# osa = HP70951B()
osa = YokogawaAQ6375()

# Initialize instruments
ps.initialize()
osa.initialize()

# Set OSA params
osa.set_rbw(RBW)
osa.set_center_wavelength(center_wav)
osa.set_ref_level(ref_level)
osa.set_span(span)

# file_name = 'LL_dev1_50V_pm1258_Solstis980nm'  # Filename where to save csv data
file_name = "waveMap_%s_%d-%dV_%dC_%s_%3.2fmW_RBW%3.2fnm" % (
    device, init_voltage, end_voltage, temp, pump_laser, pump_power*IL, RBW)  # Filename where to save csv data

# SET UP THE EXPERIMENT
instr_dict = {
    'osa': osa,
    'ps': ps
              }

params = {
    "voltage_list": voltage_list, 
    "device": device, 
    "pump_laser": pump_laser, 
    "pump_power": pump_power, 
    "temp": temp, 
    "IL": IL
    }

#%% Perform experiment
# RUN IT
[meas_volt_list, pk_wl_list] = perform_experiment(params, instrument_dict=instr_dict, filename=file_name)

# CLOSE INSTRUMENTS
osa.close()
ps.close()

#%% Plot data
fig, ax = plt.subplots(1,1, figsize=(6,4), tight_layout=True)
ax.plot(meas_volt_list, np.array(pk_wl_list)*1e9)
ax.set_xlabel("Voltage (V)")
ax.set_ylabel("Wavelength (nm)")
# %%

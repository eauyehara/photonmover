from photonmover.Interfaces.Experiment import Experiment
from photonmover.utils.plot_utils import plot_graph

# Interfaces/instruments necessary for the experiment
# - You use an Interface if any instrument of that category can be used
# - You use a specific instrument if you can only use that specific model
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

class vcsel_sweep_delay(Experiment):

    def __init__(self, instrument_list, visa_lock=None):
        """
        :param instrument_list: list of available instruments. IMPORTANT: WE ASSUME THAT THE INSTRUMENTS HAVE BEEN INITIALIZED ALREADY!
        """
        super().__init__(visa_lock)

        # It is always good practice to initialize variables in the init

        # Instruments
        self.voa = None
        self.ps = None
        self.det = None

        self.data = None

        if not self.check_necessary_instruments(instrument_list):
            raise ValueError("The necessary instruments for this experiment are not present!")

    def check_necessary_instruments(self, instrument_list):
        """
        Checks if the instruments necessary to perform the experiment are present.
        :param instrument_list: list of the available instruments
        :return: True if the necessary instruments are present, False otherwise.
        """

        for instr in instrument_list:
            if isinstance(instr, Keithley2635A):
                self.ps = instr
            if isinstance(instr, HP54750A):
                self.det = instr
            if isinstance(instr, HP70951B):
                self.det = instr
            if isinstance(instr, AgilentE3633A):
                self.voa = instr


        if (self.ps is not None) and (self.det is not None) and (self.voa is not None):
            return True
        else:
            return False

    def get_description(self):
        """
        Returns a string with a brief summary of the experiment.
        """
        return "Sweeps power supply voltage while recording osc trace into single csv file"

    def get_name(self):
        """
        Returns a string with the experiment name
        """
        return "vcsel_sweep_delay"

    def perform_experiment(self, params, filename=None):
        """
        Performs the experiment, and saves the relevant data (if there is any)
        to the specified file (if given) - assumes osa set to desired display parameters prior to running experiment
        :param params: dictionary of the parameters necessary for the experiment.
        :param filename: if specified, the data is saved in the specified file.
        :return: [meas_volt_list, pk_wl_list]
        """

        params = self.check_all_params(params)

        voltage_list = params["voltage_list"]
        voa_set = params["voa_set"]

        meas_volt_list = []
        t_peak = []

        if isinstance(self.det, HP54750A):
            # Read trace from sampling oscilloscope
            full_filename = params["detector"] + "_" + filename
            [t, _] = self.det.read_waveform([1])
            trace_len = np.array(t).shape[1]
            trace_array = np.array(np.zeros((len(voltage_list), int(trace_len))))
            x_axis = t[0]
        elif isinstance(self.det, HP70951B):
            [wavelength, _] = self.det.read_data()
            trace_len = np.array(wavelength).shape[0]
            trace_array = np.array(np.zeros((len(voltage_list), int(trace_len))))
            x_axis = wavelength

        # Sweep power supply voltage and get osa trace
        for ind, volt in enumerate(voltage_list):

            print('Setting power supply to %.4f V...' % volt)
            # Set the voltage
            self.ps.set_voltage(volt)
            self.voa.set_voltage(voa_set[ind])

            meas_volt = self.ps.measure_voltage()
            print('Power Supply voltage set to %0.4f V' % meas_volt)
            meas_volt_list.append(meas_volt) #[V]

            # Wait [s]
            # time.sleep(10)

            #Read osc trace and store in wavelength_array
            if isinstance(self.det, HP54750A):
                # Read trace from sampling oscilloscope
                full_filename = params["detector"] + "_" + filename
                [_, waveform] = self.det.read_waveform([1])
                trace_array[ind, :] = np.reshape(np.array(waveform), (1, int(trace_len)))

            elif isinstance(self.det, HP70951B):
                full_filename = params["detector"] + "_" + filename
                [_, waveform] = self.det.read_data()
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

        #Ramp back to initial voltage
        # reverse_voltage_list = list(voltage_list)
        # reverse_voltage_list.reverse()
        # reverse_voltage_str = str(reverse_voltage_list).replace('[', '{')
        # reverse_voltage_str = reverse_voltage_str.replace(']', '}')
        #
        # ps.linear_volt_sweep(volt_list=reverse_voltage_str, settling_time=.3, num_points=len(reverse_voltage_list))

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

        self.data = [meas_volt_list, t_peak]

        return [meas_volt_list, t_peak]

    def required_params(self):
        """
        Returns a list with the keys that need to be specified in the params dictionary, in order for
        a measurement to be performed
        """

        return ["voltage_list", "voa_set", "device", "pump_laser","pump_power","temp"]


    def plot_data(self, canvas_handle, data=None):

        if data is None:
            if self.data is not None:
                data = self.data
            else:
                raise ValueError('plot_data was called before performing the experiment or providing data')

        voltage = data[0]
        delay = data[1]

        plot_graph(x_data=voltage, y_data=delay, canvas_handle=canvas_handle, xlabel='Voltage (V)',
                   ylabel='Peak Delay (ps)', title='Tuning Curve', legend=None)

import sys
# from pyqtgraph.Qt import QtGui, QtCore
from PyQt5 import QtGui, QtCore, QtWidgets
import pyqtgraph as pg

# class Window(QtGui.QMainWindow):
class Window(QtWidgets.QMainWindow):
    def __init__(self):
        super(Window, self).__init__()
        self.setGeometry(100, 100, 1000, 500)
        self.setWindowTitle("Delay Sweep")

        # Menu definition
        mainMenu = self.menuBar()

        # Set Window as central widget
        # self.w = QtGui.QWidget()
        self.w = QtWidgets.QWidget()
        self.setCentralWidget(self.w)

        ## Create a grid layout to manage the widgets size and position
        # self.layout = QtGui.QGridLayout()
        self.layout = QtWidgets.QGridLayout()
        self.w.setLayout(self.layout)

        # plot widget
        self.p_power = pg.PlotWidget()
        self.xlabel = self.p_power.setLabel('bottom', text='Voltage', units='V')
        self.ylabel = self.p_power.setLabel('left', text='Delay', units='ps')
        self.layout.addWidget(self.p_power, 0, 0)

        self.p_power_handle = self.p_power.plot(pen=(1, 3))

        self.show()


    def plot(self, x, y):
        self.p_power_handle.setData(x,y)

if __name__ == '__main__':
    # ------------------------------------------------------------
    # SAFETY LIMITS
    i_limit = 100e-9  # current limit

    # OTHER PARAMETERS
    detector = 'osa'
    device = 'dev1b_delayvolt13'
    pump_laser = 'OE1076CW115.7mA' #'OEland1076' #'OEland1038' #'CW976'
    pump_power = 0 #mW
    IL = 0.62
    RBW = 0.1 #nm
    temp = 15 #C

    # EXPERIMENT PARAMETERS
    init_voltage = 0  # [V]
    end_voltage = 72 # [V]
    increment = 1  # Voltage increment
    voltage_list = np.arange(init_voltage, end_voltage+increment, increment) #end_voltage+1 or will stop at end_voltage-1
    # voltage_list = np.append(voltage_list, [72.2])
    print(voltage_list)
    delayvolt_file = "delayvolt13"
    delay_dir = os.path.join(os.getcwd(),"delayvolt")
    # ------------------------------------------------------------

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
    instr_list = [det, ps, voa]
    params = {"detector": detector, "voltage_list": voltage_list, "voa_set": voa_set, "device": device, "pump_laser": pump_laser, "pump_power": pump_power, "temp": temp}
    exp = vcsel_sweep_delay(instr_list)

    # RUN IT
    [meas_volt_list, t_peak] = exp.perform_experiment(params, filename=file_name)

    # PLOT DATA
    # app = QtGui.QApplication(sys.argv)
    app = QtWidgets.QApplication(sys.argv)
    GUI = Window()
    GUI.plot(meas_volt_list, t_peak)

    # CLOSE INSTRUMENTS
    det.close()
    ps.close()

    sys.exit(app.exec_())
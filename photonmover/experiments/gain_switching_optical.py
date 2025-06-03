from photonmover.Interfaces.Experiment import Experiment
from photonmover.utils.plot_utils import plot_graph

# Interfaces/instruments necessary for the experiment
# - You use an Interface if any instrument of that category can be used
# - You use a specific instrument if you can only use that specific model
from photonmover.instruments.Power_meters.Thorlabs import ThorlabsPowerMeter
from photonmover.instruments.Power_Supplies.AgilentE3633A import AgilentE3633A
from photonmover.instruments.Oscilloscopes.HP54750A import HP54750A
from photonmover.instruments.Optical_spectrum_analyzers.HP70951B import HP70951B

# General imports
import time
import numpy as np
import csv
from scipy import io

from ctypes import create_string_buffer, c_int16, c_double, byref, c_char_p, c_bool


class gain_switching_optical(Experiment):

    def __init__(self, instrument_list, visa_lock=None):
        """
        :param instrument_list: list of available instruments. IMPORTANT: WE ASSUME THAT THE INSTRUMENTS HAVE BEEN INITIALIZED ALREADY!
        """
        super().__init__(visa_lock)

        # It is always good practice to initialize variables in the init

        # Instruments
        self.pm = None
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
            if isinstance(instr, ThorlabsPowerMeter):
                self.pm = instr
            if isinstance(instr, AgilentE3633A):
                self.ps = instr
            if isinstance(instr, HP54750A):
                self.det = instr
            if isinstance(instr, HP70951B):
                self.det = instr

        if (self.pm is not None) and (self.ps is not None) and (self.det is not None):
            return True
        else:
            return False

    def get_description(self):
        """
        Returns a string with a brief summary of the experiment.
        """
        return "Gain switching experiment for Praevium optically pumped VCSELs - ramps VOA bias to vary optical pumping power."

    def get_name(self):
        """
        Returns a string with the experiment name
        """
        return "gain switching optical"
        
    def perform_experiment(self, params, filename=None):
        """
        Performs the experiment, and saves the relevant data (if there is any)
        to the specified file (if given)
        :param params: dictionary of the parameters necessary for the experiment.
        :param filename: if specified, the data is saved in the specified file.
        :return:
        """

        params = self.check_all_params(params)

        fiber_coupling_efficiency = params["fiber_coupling_efficiency"]  #collection bench fiber coupling
        pump_wavelength = params["pump_wavelength"] #pump wavelength [nm]
        num_avg = params["num_avg"]  # sampling oscilloscope
        tuning_voltage = params["tuning_voltage"]  #VCSEL tuning voltage
        VOA_voltage = params["VOA_voltage"] #VOA voltage [V]
        IL = params["IL"] #Insertion loss of WDM and fiber connection to 99:1 tap
        RBW = params["RBW"] #OSA RBW

        # Set power meter wavelength
        attribute = c_int16(0) #set wavelength
        meas_pump_wavelength = c_double()

        self.pm.getWavelength(attribute, byref(meas_pump_wavelength))
        if meas_pump_wavelength is not pump_wavelength:
            self.pm.setWavelength(c_double(pump_wavelength))

        pump_power_list = []
        trace_peak_list = []
        meas_volt_list = []

        # Sweep VOA voltage and get power
        for volt in VOA_voltage:

            print('Setting VOA to %.4f V...' % volt)
            # Set the voltage
            self.ps.set_voltage(volt)

            meas_volt = self.ps.measure_voltage()
            print('VOA voltage set to %0.4f V' % meas_volt)
            meas_volt_list.append(meas_volt) #[V]

            meas_curr = self.ps.measure_current()
            print('VOA current is %0.4f A' % meas_curr)

            # Wait time (increase if oscilloscope num_avg or number of points large)
            time.sleep(10.0)  # [s]

            #Read pump power through 1% tap, scale to power in 99% tap
            [pump_power_tap, _] = self.pm.get_powers()
            pump_power = pump_power_tap*99*IL  #scale to full power
            pump_power_list.append(pump_power)  # [W]
            print('Pump power at %0.6f mW' % (pump_power*1e3))

            # Create full filename for oscilloscope trace
            if filename is not None:
                # Save the data in a csv file
                time_tuple = time.localtime()
                full_filename = "%s_%.1fmW_%.1fV_%d-%d-%d" % (
                    filename, pump_power * 1e3, volt,
                    time_tuple[0],
                    time_tuple[1],
                    time_tuple[2])
                print(full_filename)

            if isinstance(self.det, HP54750A):
                # Read trace from sampling oscilloscope
                det = "osc"
                full_filename = det + "_" + full_filename
                [_, waveform] = self.det.read_waveform([1], full_filename)

            elif isinstance(self.det, HP70951B):
                det = "osa"
                full_filename = det + "_" + full_filename
                [_, waveform] = self.det.read_data(filename=full_filename)


            trace_peak_list.append(max(waveform))
            #print('Trace peak voltage is %0.6f V' % (max(waveform)*1e3))


        print('Finished gain switching experiment')
        print('-----------------------------')

        # Save the parameters in a .mat file
        time_tuple = time.localtime()
        params_filename = "%s_params_%s_%.1f-%0.1fmW_%d-%d-%d.mat" % (
        det, filename, pump_power_list[0] * 1e3, pump_power_list[-1] * 1e3,
        time_tuple[0],
        time_tuple[1],
        time_tuple[2])

        io.savemat(params_filename, {'pump_power_list': pump_power_list,
                                     'VOA_voltage_list': meas_volt_list,
                                     'pump_wavelength': pump_wavelength,
                                     'tuning_voltage': tuning_voltage,
                                     'fiber_coupling_efficiency': fiber_coupling_efficiency,
                                     'IL': IL,
                                     'RBW': RBW,
                                     'num_avg': num_avg})
        # units: pump_power_list: [W], VOA_voltage_list: [V], pump_wavelength: [nm], tuning_voltage: [V], fiber_coupling efficiency: [], 'num_avg': [#]

        self.data = [pump_power_list, trace_peak_list]
        print(np.array(pump_power_list).shape)
        print(np.array(np.array(trace_peak_list)))
        return [np.array(pump_power_list), np.array(trace_peak_list)]
    
    def required_params(self):
        """
        Returns a list with the keys that need to be specified in the params dictionary, in order for
        a measurement to be performed
        """
        return ["fiber_coupling_efficiency","IL", "RBW", "pump_wavelength", "num_avg", "tuning_voltage", "VOA_voltage"]

    def plot_data(self, canvas_handle, data=None):
        
        if data is None:
            if self.data is not None:
                data = self.data
            else:
                raise ValueError('plot_data was called before performing the experiment or providing data')
        
        pump_power = data[0]
        trace_peak = data[1]
        if isinstance(self.det, HP54750A):
            plot_graph(x_data=pump_power, y_data=trace_peak, canvas_handle=canvas_handle, xlabel='Pump Power (mW)', ylabel='Trace Peak (V)', title='Peak Trace Voltage', legend=None)
        elif isinstance(self.det, HP70951B):
            plot_graph(x_data=pump_power, y_data=trace_peak, canvas_handle=canvas_handle, xlabel='Pump Power (mW)', ylabel='Spectrum Peak (dBm)', title='Peak Amplitiude', legend=None)


import sys
# from pyqtgraph.Qt import QtGui, QtCore
from PyQt5 import QtGui, QtCore, QtWidgets
import pyqtgraph as pg

# class Window(QtGui.QMainWindow):
class Window(QtWidgets.QMainWindow):
    def __init__(self):
        super(Window, self).__init__()
        self.setGeometry(100, 100, 1000, 500)
        self.setWindowTitle("PER Sweep")

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
        self.xlabel = self.p_power.setLabel('bottom', text='Pump Power', units='W')
        self.ylabel = self.p_power.setLabel('left', text='Peak Amplitude', units='a.u.')
        self.layout.addWidget(self.p_power, 0, 0)

        self.p_power_handle = self.p_power.plot(pen=(1, 3))

        self.show()


    def plot(self, x, y):
        self.p_power_handle.setData(x,y)

# import sys
# from pyqtgraph.Qt import QtGui, QtCore
# import pyqtgraph as pg
#
# class Window(QtGui.QMainWindow):
#     def __init__(self):
#         super(Window, self).__init__()
#         self.setGeometry(100, 100, 1000, 500)
#         self.setWindowTitle("Peak Trace Voltage")
#
#         # Menu definition
#         mainMenu = self.menuBar()
#
#         # Set Window as central widget
#         self.w = QtGui.QWidget()
#         self.setCentralWidget(self.w)
#
#         ## Create a grid layout to manage the widgets size and position
#         self.layout = QtGui.QGridLayout()
#         self.w.setLayout(self.layout)
#
#         # plot widget
#         self.p_power = pg.PlotWidget()
#         self.xlabel = self.p_power.setLabel('bottom', text='Pump Power', units='W')
#         self.ylabel = self.p_power.setLabel('left', text='Peak Amplitude', units='a.u.')
#         self.layout.addWidget(self.p_power, 0, 0)
#
#         self.p_power_handle = self.p_power.plot(pen=(1, 3))
#
#         self.show()
#
#
#     def plot(self, x, y):
#         self.p_power_handle.setData(x,y)

if __name__ == '__main__':

    # -------------------------------------------------------------
    # SAFETY LIMITS
    i_limit = 10e-9#current limit

    # POWER METER SETTINGS
    pump_wavelength = 1076  #nm

    # SPECIFY DETECTOR ("osa" or "osc")
    detector = "osc"

    # OSCILLOSCOPE PARAMETERS
    num_avg = 16

    # OSA PARAMETERS
    RBW = 0.1 #nm

    #DEVICE PARAMETERS
    device = 'dev1b-CW92mA'
    tuning_voltage = 88.6   # [V]

    # COLLECTION BENCH PARAMETERS
    fiber_coupling_efficiency = 1
    IL = 0.54 #insertion loss measured as (WDM 980/1310 output) / (1% tap)*100 - scales 1% tap output to actual input to VCSEL fiber

    # EXPERIMENT PARAMETERS
    init_voltage = 3.1 #2.5 #2.8 #2.46 #.99  # [V] Minimum transmission on VOA (Note: when set to 5V, AgilentE3633A momentarily exceeds current limit when turning output on)
    end_voltage = 2.4 #1.8 #2.2 #1.6  # [V] Maximum transmission on VOA
    num_points = 15  # Number of points between init and end current
    VOA_voltage = np.linspace(init_voltage, end_voltage, num_points)

    # -----------------------------------------------------------

    # INSTRUMENTS
    ps = AgilentE3633A(current_limit=i_limit)
    pm = ThorlabsPowerMeter()
    if detector == "osc":
        det = HP54750A()
    elif detector == "osa":
        det = HP70951B()


    # Set up power supply
    ps.initialize()

    # Set up Thorlabs power meter
    pm.initialize()
    pm.setPowerAutoRange(c_int16(1))  #enable autorange

    # Set up detector
    det.initialize()

    #file_name = 'dev1_50V_1064nm'  # Filename where to save csv data
    file_name = "%s_%dV_%dnm" % (device, tuning_voltage, pump_wavelength) # Filename where to save csv data

    # SET UP THE EXPERIMENT
    instr_list = [pm, ps, det]

    params = {"fiber_coupling_efficiency": fiber_coupling_efficiency, "IL": IL, "RBW": RBW, "pump_wavelength": pump_wavelength, "num_avg": num_avg, "tuning_voltage": tuning_voltage, "VOA_voltage": VOA_voltage}
    exp = gain_switching_optical(instr_list)

    # RUN IT
    [pump_power_list, trace_peak_list] = exp.perform_experiment(params, filename=file_name)


    # PLOT DATA
    # app = QtGui.QApplication(sys.argv)
    app = QtWidgets.QApplication(sys.argv)
    GUI = Window()
    GUI.plot(pump_power_list, trace_peak_list)
    sys.exit(app.exec_())

    # CLOSE INSTRUMENTS
    pm.close()
    ps.close()
    det.close()

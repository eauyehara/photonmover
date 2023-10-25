from photonmover.Interfaces.Experiment import Experiment
from photonmover.utils.plot_utils import plot_graph

# Interfaces/instruments necessary for the experiment
# - You use an Interface if any instrument of that category can be used
# - You use a specific instrument if you can only use that specific model
from photonmover.instruments.Power_meters.Thorlabs import ThorlabsPowerMeter
from photonmover.instruments.Oscilloscopes.HP54750A import HP54750A
from photonmover.instruments.Optical_spectrum_analyzers.HP70951B import HP70951B
from photonmover.instruments.Laser_drivers.ITC4001 import ITC4001

# General imports
import time
import numpy as np
import csv
from scipy import io

from ctypes import create_string_buffer, c_int16, c_double, byref, c_char_p, c_bool


class BOA_current_sweep(Experiment):

    def __init__(self, instrument_list, visa_lock=None):
        """
        :param instrument_list: list of available instruments. IMPORTANT: WE ASSUME THAT THE INSTRUMENTS HAVE BEEN INITIALIZED ALREADY!
        """
        super().__init__(visa_lock)

        # It is always good practice to initialize variables in the init

        # Instruments
        self.pm = None
        self.ld = None
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
            if isinstance(instr, ITC4001):
                self.ld = instr
            if isinstance(instr, HP54750A):
                self.det = instr
            if isinstance(instr, HP70951B):
                self.det = instr

        if (self.pm is not None) and (self.ld is not None) and (self.det is not None):
            return True
        else:
            return False

    def get_description(self):
        """
        Returns a string with a brief summary of the experiment.
        """
        return "Measure sampling scope (or OSA) trace while sweeping BOA current"

    def get_name(self):
        """
        Returns a string with the experiment name
        """
        return "BOA current sweep"
        
    def perform_experiment(self, params, filename=None):
        """
        Performs the experiment, and saves the relevant data (if there is any)
        to the specified file (if given)
        :param params: dictionary of the parameters necessary for the experiment.
        :param filename: if specified, the data is saved in the specified file.
        :return:
        """

        params = self.check_all_params(params)

        pump_wavelength = params["pump_wavelength"] #pump wavelength [nm]
        num_avg = params["num_avg"]  # sampling oscilloscope
        tuning_voltage = params["tuning_voltage"]  #VCSEL tuning voltage
        IL = params["IL"] #Insertion loss of WDM and fiber connection to 99:1 tap
        RBW = params["RBW"] #OSA RBW
        BOA_curr = params["BOA_curr"] #BOA current

        # Set power meter wavelength
        attribute = c_int16(0) #set wavelength
        meas_pump_wavelength = c_double()

        self.pm.getWavelength(attribute, byref(meas_pump_wavelength))
        if meas_pump_wavelength is not pump_wavelength:
            self.pm.setWavelength(c_double(pump_wavelength))

        trace_peak_list = []
        meas_curr_list = []

        # Read pump power through 1% tap, scale to power in 99% tap
        [pump_power_tap, _] = self.pm.get_powers()
        pump_power = pump_power_tap * 100 * IL  # scale to full power
        print('Pump power at %0.6f mW' % (pump_power * 1e3))

        # Sweep VOA voltage and get power
        for curr in BOA_curr:
            print('Setting BOA current to %.4f A...' % curr)
            # Set the BOA current
            self.ld.set_current(curr)

            # Measure BOA current
            time.sleep(0.5) #[s]
            meas_curr = self.ld.measure_current()
            print('BOA current set to %0.4f A' % meas_curr)
            meas_curr_list.append(meas_curr) #[A]

            # Wait time (increase if oscilloscope num_avg large)
            time.sleep(7.0)  # [s]

            # Create full filename for oscilloscope trace
            if filename is not None:
                # Save the data in a csv file
                time_tuple = time.localtime()
                full_filename = "BOAsweep_%s_%.1fmW_%.1fmA_%d-%d-%d" % (
                    filename, pump_power * 1e3, curr*1e3,
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
        params_filename = "params_%s_BOAsweep_%s_%.1f-%0.1fmA_%d-%d-%d.mat" % (
        det, filename, meas_curr_list[0] * 1e3, meas_curr_list[-1] * 1e3,
        time_tuple[0],
        time_tuple[1],
        time_tuple[2])

        io.savemat(params_filename, {'pump_power': pump_power,
                                     'meas_curr_list': meas_curr_list,
                                     'pump_wavelength': pump_wavelength,
                                     'tuning_voltage': tuning_voltage,
                                     'VCSEL_wavelength': VCSEL_wavelength,
                                     'IL': IL,
                                     'RBW': RBW,
                                     'num_avg': num_avg})
        # units: pump_power: [W], meas_curr_list: [A], pump_wavelength: [nm], tuning_voltage: [V], IL [], RBW [nm], num_avg: [#]

        self.data = [np.array(meas_curr_list), np.array(trace_peak_list)]
        return [np.array(meas_curr_list), np.array(trace_peak_list)]
    
    def required_params(self):
        """
        Returns a list with the keys that need to be specified in the params dictionary, in order for
        a measurement to be performed
        """
        return ["IL", "RBW", "pump_wavelength", "VCSEL_wavelength", "num_avg", "tuning_voltage", "BOA_curr"]

    def plot_data(self, canvas_handle, data=None):
        
        if data is None:
            if self.data is not None:
                data = self.data
            else:
                raise ValueError('plot_data was called before performing the experiment or providing data')
        
        BOA_current = data[0]*1e3
        trace_peak = data[1]
        if isinstance(self.det, HP54750A):
            plot_graph(x_data=BOA_current, y_data=trace_peak, canvas_handle=canvas_handle, xlabel='BOA Current (mA)', ylabel='Trace Peak (V)', title='Peak Trace Voltage', legend=None)
        elif isinstance(self.det, HP70951B):
            plot_graph(x_data=BOA_current, y_data=trace_peak, canvas_handle=canvas_handle, xlabel='BOA Current (mA)', ylabel='Spectrum Peak (dBm)', title='Peak Amplitiude', legend=None)




import sys
from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg

class Window(QtGui.QMainWindow):
    def __init__(self):
        super(Window, self).__init__()
        self.setGeometry(100, 100, 1000, 500)
        self.setWindowTitle("Peak Trace Voltage")

        # Menu definition
        mainMenu = self.menuBar()

        # Set Window as central widget
        self.w = QtGui.QWidget()
        self.setCentralWidget(self.w)

        ## Create a grid layout to manage the widgets size and position
        self.layout = QtGui.QGridLayout()
        self.w.setLayout(self.layout)

        # plot widget
        self.p_power = pg.PlotWidget()
        self.xlabel = self.p_power.setLabel('bottom', text='BOA Current', units='mA')
        self.ylabel = self.p_power.setLabel('left', text='Peak Amplitude', units='a.u.')
        self.layout.addWidget(self.p_power, 0, 0)

        self.p_power_handle = self.p_power.plot(pen=(1, 3))

        self.show()


    def plot(self, x, y):
        self.p_power_handle.setData(x,y)

if __name__ == '__main__':

    # -------------------------------------------------------------
    # SAFETY LIMITS
    i_limit = 0.705 #[A] current limit

    # POWER METER SETTINGS
    pump_wavelength = 1040  #nm

    # SPECIFY DETECTOR ("osa" or "osc")
    detector = "osa"

    # OSCILLOSCOPE PARAMETERS
    num_avg = 64

    # OSA PARAMETERS
    RBW = 1 #nm

    #DEVICE PARAMETERS
    device = 'bf_PolStable'
    tuning_voltage = 45  # [V]
    VCSEL_wavelength = 1250 #[nm]

    # COLLECTION BENCH PARAMETERS
    IL = 0.69 #insertion loss measured as (WDM 980/1310 output) / (1% tap)*100 - scales 1% tap output to actual input to VCSEL fiber

    # EXPERIMENT PARAMETERS
    BOA_curr = np.concatenate(([150], np.arange(180,260,10), np.arange(300,750,50)))*1e-3 #[A]
    # BOA_curr = np.array([0.01, 0.02, 0.03, 0.04, 0.05])

    # ------------------------------------------------------------

    # INSTRUMENTS
    ld = ITC4001()
    pm = ThorlabsPowerMeter()
    if detector == "osc":
        det = HP54750A()
    elif detector == "osa":
        det = HP70951B()


    # Set up laser driver (for BOA)
    ld.initialize()
    ld.set_current_limit(current_limit=i_limit)

    # Set up Thorlabs power meter
    pm.initialize()
    pm.setPowerAutoRange(c_int16(1))  #enable autorange

    # Set up detector
    det.initialize()

    file_name = "%s_%dnm_%dnm" % (device, VCSEL_wavelength, pump_wavelength) # Filename where to save csv data

    # SET UP THE EXPERIMENT
    instr_list = [pm, ld, det]

    params = {"IL": IL, "RBW": RBW, "pump_wavelength": pump_wavelength, "num_avg": num_avg, "tuning_voltage": tuning_voltage, "VCSEL_wavelength": VCSEL_wavelength, "BOA_curr": BOA_curr}
    exp = BOA_current_sweep(instr_list)

    # RUN IT
    [curr_list, trace_peak_list] = exp.perform_experiment(params, filename=file_name)


    # PLOT DATA
    app = QtGui.QApplication(sys.argv)
    GUI = Window()
    GUI.plot(curr_list, trace_peak_list)
    sys.exit(app.exec_())

    # CLOSE INSTRUMENTS
    pm.close()
    ld.close()
    det.close()

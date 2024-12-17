
from photonmover.Interfaces.Experiment import Experiment
from photonmover.utils.plot_utils import plot_graph

# Interfaces/instruments necessary for the experiment
# - You use an Interface if any instrument of that category can be used
# - You use a specific instrument if you can only use that specific model
from photonmover.instruments.Power_meters.Thorlabs import ThorlabsPowerMeter
from photonmover.instruments.Power_Supplies.AgilentE3633A import AgilentE3633A
from photonmover.instruments.Multimeters.HP34401 import HP34401

# General imports
import time
import numpy as np
import csv
from scipy import io

from ctypes import create_string_buffer, c_int16, c_double, byref, c_char_p, c_bool


class PD_saturation(Experiment):

    def __init__(self, instrument_list, visa_lock=None):
        """
        :param instrument_list: list of available instruments. IMPORTANT: WE ASSUME THAT THE INSTRUMENTS HAVE BEEN INITIALIZED ALREADY!
        """
        super().__init__(visa_lock)

        # It is always good practice to initialize variables in the init

        # Instruments
        self.ps = None
        self.pm = None
        self.mm = None

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
            if isinstance(instr, AgilentE3633A):
                self.ps = instr
            if isinstance(instr, ThorlabsPowerMeter):
                self.pm = instr
            if isinstance(instr, HP34401):
                self.mm = instr

        if (self.ps is not None) and (self.pm is not None) and (self.mm is not None):
            return True
        else:
            return False

    def get_description(self):
        """
        Returns a string with a brief summary of the experiment.
        """
        return "Measures power saturation of photodiode by ramping input power while monitoring voltage across resistor"

    def get_name(self):
        """
        Returns a string with the experiment name
        """
        return "PD_saturation"

    def perform_experiment(self, params, filename=None):
        """
        Performs the experiment, and saves the relevant data (if there is any)
        to the specified file (if given) - assumes osa set to desired display parameters prior to running experiment
        :param params: dictionary of the parameters necessary for the experiment.
        :param filename: if specified, the data is saved in the specified file.
        :return: [voltage_list, peak_wl_list, peak_dBm_list, power_list]
        """

        params = self.check_all_params(params)

        voltage_list = params["voltage_list"]
        IL = params["IL"]
        pump_wavelength = params["pump_wavelength"]
        RL = params["load_resistor"]
        Rf = params["filter_resistor"]
        Cf = params["filter_capacitor"]

        VOA_volt_list = []
        power_list = []
        PD_volt_list = []


        # Sweep power supply voltage to VOA, read power meter tap, read multimeter voltage
        for ind, volt in enumerate(voltage_list):

            print('Setting power supply to %.4f V...' % volt)
            # Set the VOA voltage
            self.ps.set_voltage(volt)

            # Read actual VOA voltage
            meas_volt = self.ps.measure_voltage()
            print('Power Supply voltage set to %0.4f V' % meas_volt)
            VOA_volt_list.append(meas_volt) #[V]

            # Wait [s]
            time.sleep(0.5)

            # Read optical power through tap
            [tap_power, _] = self.pm.get_powers()
            power_list.append(tap_power*100*IL)
            # power_list.append(tap_power*2)

            # Read PD voltage from multimeter
            PD_volt = self.mm.measure_voltage('DC')
            PD_volt_list.append(PD_volt)

        print('Finished voltage sweep')
        print('-----------------------------')
        self.mm.beep()

        if filename is not None:
            # Save the data in a csv file
            time_tuple = time.localtime()
            complete_filename = "./data/%s_%d-%d-%d_%d-%d-%d.csv" % (filename,
                                                            time_tuple[0],
                                                            time_tuple[1],
                                                            time_tuple[2],
                                                            time_tuple[3],
                                                            time_tuple[4],
                                                            time_tuple[5])

            with open(complete_filename, 'w+') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(VOA_volt_list)  #[V]
                writer.writerow(power_list)  #[W]
                writer.writerow(PD_volt_list) #[V]

                # Save the parameters in a .mat file
                time_tuple = time.localtime()
                params_filename = "./data/params_%s_%d-%d-%d_%d-%d-%d.mat" % (
                    filename,
                    time_tuple[0],
                    time_tuple[1],
                    time_tuple[2],
                    time_tuple[3],
                    time_tuple[4],
                    time_tuple[5])

                io.savemat(params_filename, {'voltage_list': voltage_list,
                                             'IL': IL,
                                             'pump_wavelength': pump_wavelength,
                                             'load_resistor': RL,
                                             'filter_resitor': Rf,
                                             'filter_capacitor': Cf
                                             })

        self.data = [power_list, PD_volt_list]

        return [power_list, PD_volt_list]

    def required_params(self):
        """
        Returns a list with the keys that need to be specified in the params dictionary, in order for
        a measurement to be performed
        """
        return ["voltage_list", "IL", "pump_wavelength", "load_resistor", "filter_resistor", "filter_capacitor"]

    def plot_data(self, canvas_handle, data=None):

        if data is None:
            if self.data is not None:
                data = self.data
            else:
                raise ValueError('plot_data was called before performing the experiment or providing data')

        power_list = data[0]
        PD_voltage = data[1]

        plot_graph(x_data=power_list, y_data=PD_voltage, canvas_handle=canvas_handle, xlabel='Power (W)',
                   ylabel='Voltage (V)', title='PD Saturation', legend=None)

import sys
# from pyqtgraph.Qt import QtGui, QtCore
from PyQt5 import QtGui, QtCore, QtWidgets
import pyqtgraph as pg

# class Window(QtGui.QMainWindow):
class Window(QtWidgets.QMainWindow):
    def __init__(self):
        super(Window, self).__init__()
        self.setGeometry(100, 100, 1000, 500)
        self.setWindowTitle("PD Saturation")

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
        self.xlabel = self.p_power.setLabel('bottom', text='Power', units='W')
        self.ylabel = self.p_power.setLabel('left', text='Photodiode Voltage', units='V')
        self.layout.addWidget(self.p_power, 0, 0)

        self.p_power_handle = self.p_power.plot(pen=(1, 3))

        self.show()


    def plot(self, x, y):
        self.p_power_handle.setData(x,y)



if __name__ == '__main__':
    # ------------------------------------------------------------
    # SAFETY LIMITS
    i_limit = 0.003  # current limit [A]

    # OTHER PARAMETERS
    photodiode = 'FGA21'
    pump_wavelength = 1038#[nm]
    IL = 0.671
    reverse_voltage = 55 #[V]
    RL = 50 #[ohms] Load resistor
    Rf = 1000 #[ohms] Filter resistor
    Cf = 0.1e-6 #[F] Filter capacitance

    # EXPERIMENT PARAMETERS
    init_voltage = 4.0  # [V] Minimum transmission on VOA (Note: when set to 5V, AgilentE3633A momentarily exceeds current limit when turning output on)
    end_voltage = 1.0  # [V] Maximum transmission on VOA
    num_points = 100  # Number of points between init and end current
    voltage_list = np.linspace(init_voltage, end_voltage, num_points)
    # ------------------------------------------------------------

    # INSTRUMENTS
    ps = AgilentE3633A(current_limit=i_limit)
    pm = ThorlabsPowerMeter()
    mm = HP34401()

    # Initialize instruments
    ps.initialize()
    pm.initialize()
    pm.setPowerAutoRange(c_int16(1))  #enable autorange
    pm.setWavelength(c_double(pump_wavelength))
    mm.initialize()


    file_name = "PDsat_%s_%dV_%dnm_%dohms" % (
        photodiode, reverse_voltage, pump_wavelength, RL)  # Filename where to save csv data

    # SET UP THE EXPERIMENT
    instr_list = [ps, pm, mm]
    params = {"voltage_list": voltage_list, "IL": IL, "pump_wavelength": pump_wavelength, "load_resistor": RL, "filter_resistor": Rf, "filter_capacitor": Cf}
    exp = PD_saturation(instr_list)

    # RUN IT
    [power_list, PD_voltage] = exp.perform_experiment(params, filename= file_name)

    # CLOSE INSTRUMENTS
    mm.close()
    ps.close()
    pm.close()

    # PLOT DATA
    # app = QtGui.QApplication(sys.argv)
    app = QtWidgets.QApplication(sys.argv)
    GUI = Window()
    GUI.plot(power_list, PD_voltage)
    sys.exit(app.exec_())
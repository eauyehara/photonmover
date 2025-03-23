from photonmover.Interfaces.Experiment import Experiment
from photonmover.utils.plot_utils import plot_graph

# Interfaces/instruments necessary for the experiment
# - You use an Interface if any instrument of that category can be used
# - You use a specific instrument if you can only use that specific model
from photonmover.instruments.Optical_spectrum_analyzers.HP70951B import HP70951B
# from photonmover.instruments.Power_Supplies.KeysightE36106B import KeysightE36106B
from photonmover.instruments.Source_meters.Keithley2635A import Keithley2635A

# General imports
import time
import numpy as np
import csv
from scipy import io

class praevium_map_wavelength(Experiment):

    def __init__(self, instrument_list, visa_lock=None):
        """
        :param instrument_list: list of available instruments. IMPORTANT: WE ASSUME THAT THE INSTRUMENTS HAVE BEEN INITIALIZED ALREADY!
        """
        super().__init__(visa_lock)

        # It is always good practice to initialize variables in the init

        # Instruments
        self.osa = None
        self.ps = None

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
            if isinstance(instr, HP70951B):
                self.osa = instr
            # if isinstance(instr, KeysightE36106B):
            if isinstance(instr, Keithley2635A):
                self.ps = instr

        if (self.osa is not None) and (self.ps is not None):
            return True
        else:
            return False

    def get_description(self):
        """
        Returns a string with a brief summary of the experiment.
        """
        return "Sweeps power supply voltage while recording osa trace into single csv file"

    def get_name(self):
        """
        Returns a string with the experiment name
        """
        return "praevium_map_wavelength"

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

        meas_volt_list = []
        pk_wl_list = []
        [trace_len, osa_settings_list] = self.osa.get_osa_parameters()

        spectrum_array = np.array(np.zeros((len(voltage_list), int(trace_len))))
        [wavelength_list, _] = self.osa.read_data()

        # Sweep power supply voltage and get osa trace
        for ind, volt in enumerate(voltage_list):

            print('Setting power supply to %.4f V...' % volt)
            # Set the voltage
            self.ps.set_voltage(volt)

            meas_volt = self.ps.measure_voltage()
            print('Power Supply voltage set to %0.4f V' % meas_volt)
            meas_volt_list.append(meas_volt) #[V]

            # Wait [s]
            time.sleep(0.5)

            #Read osa spectrum and store in wavelength_array
            [_, amp] = self.osa.read_data()
            spectrum_array[ind, :] = np.reshape(np.array(amp), (1, int(trace_len)))
            pk_ind = np.argmax(spectrum_array[ind,:])
            pk_wl_list.append(wavelength_list[pk_ind])

        print(np.shape(spectrum_array))
        print('Finished voltage sweep')
        print('-----------------------------')

        #Ramp back to initial voltage
        reverse_voltage_list = list(voltage_list)
        reverse_voltage_list.reverse()
        reverse_voltage_str = str(reverse_voltage_list).replace('[', '{')
        reverse_voltage_str = reverse_voltage_str.replace(']', '}')

        ps.linear_volt_sweep(volt_list=reverse_voltage_str, settling_time=.3, num_points=len(reverse_voltage_list))


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
                                             'voltage_list': voltage_list
                                             })

                wavvolt_filename = "%s_%dC_%s_%3.2fmW" % (
                    params["device"],
                    params["temp"],
                    params["pump_laser"],
                    params["pump_power"]*params["IL"]
                )

                self.save_wavvolt(wavvolt_filename, np.array(meas_volt_list), np.array(pk_wl_list))

        self.data = [meas_volt_list, pk_wl_list]

        return [meas_volt_list, pk_wl_list]

    def required_params(self):
        """
        Returns a list with the keys that need to be specified in the params dictionary, in order for
        a measurement to be performed
        """

        return ["voltage_list", "device", "pump_laser","pump_power","temp", "IL"]

    def save_wavvolt(self, filename, volt, peak_wl, f_interp=5000):
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


    def plot_data(self, canvas_handle, data=None):

        if data is None:
            if self.data is not None:
                data = self.data
            else:
                raise ValueError('plot_data was called before performing the experiment or providing data')

        voltage = data[0]
        peak_wl = data[1]

        plot_graph(x_data=voltage, y_data=peak_wl, canvas_handle=canvas_handle, xlabel='Voltage (V)',
                   ylabel='Peak Wavelength (nm)', title='Tuning Curve', legend=None)

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
        self.xlabel = self.p_power.setLabel('bottom', text='Voltage', units='V')
        self.ylabel = self.p_power.setLabel('left', text='Wavelength', units='nm')
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
    device = 'Dev1a'
    pump_laser = 'OEland1076' #'OEland1038' #'CW976'
    pump_power = 11 #mW
    IL = 0.52
    RBW = 0.1 #nm
    temp = 25 #C

    # EXPERIMENT PARAMETERS
    init_voltage = 0  # [V]
    end_voltage = 70 # [V]
    increment = 1  # Voltage increment
    voltage_list = np.arange(init_voltage, end_voltage+increment, increment) #end_voltage+1 or will stop at end_voltage-1
    # ------------------------------------------------------------

    # INSTRUMENTS
    # ps = KeysightE36106B(current_limit=i_limit)
    ps = Keithley2635A(current_compliance=100e-9, voltage_compliance=81) #A, V
    osa = HP70951B()

    # Initialize instruments
    ps.initialize()
    osa.initialize()

    # file_name = 'LL_dev1_50V_pm1258_Solstis980nm'  # Filename where to save csv data
    file_name = "waveMap_%s_%d-%dV_%dC_%s_%3.2fmW_RBW%3.2fnm" % (
        device, init_voltage, end_voltage, temp, pump_laser, pump_power*IL, RBW)  # Filename where to save csv data

    # SET UP THE EXPERIMENT
    instr_list = [osa, ps]
    params = {"voltage_list": voltage_list, "device": device, "pump_laser": pump_laser, "pump_power": pump_power, "temp": temp, "IL": IL}
    exp = praevium_map_wavelength(instr_list)

    # RUN IT
    [meas_volt_list, pk_wl_list] = exp.perform_experiment(params, filename=file_name)

    # PLOT DATA
    # app = QtGui.QApplication(sys.argv)
    app = QtWidgets.QApplication(sys.argv)
    GUI = Window()
    GUI.plot(meas_volt_list, pk_wl_list)

    # CLOSE INSTRUMENTS
    osa.close()
    ps.close()

    sys.exit(app.exec_())
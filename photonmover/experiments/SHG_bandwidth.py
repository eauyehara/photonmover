from photonmover.Interfaces.Experiment import Experiment
from photonmover.utils.plot_utils import plot_graph

# Interfaces/instruments necessary for the experiment
# - You use an Interface if any instrument of that category can be used
# - You use a specific instrument if you can only use that specific model
from photonmover.instruments.Optical_spectrum_analyzers.HP70951B import HP70951B
from photonmover.instruments.Source_meters.Keithley2635A import Keithley2635A
from photonmover.instruments.Power_meters.Thorlabs import ThorlabsPowerMeter

# General imports
import time
import numpy as np
import csv
from scipy import io
from ctypes import create_string_buffer, c_int16, c_double, byref, c_char_p


#Compatible Power meter sensors
SENSOR_1 = b'S150C'  # 350-1100nm
SENSOR_2 = b'S130C'  # 400-1100nm

class SHG_bandwidth(Experiment):

    def __init__(self, instrument_list, visa_lock=None):
        """
        :param instrument_list: list of available instruments. IMPORTANT: WE ASSUME THAT THE INSTRUMENTS HAVE BEEN INITIALIZED ALREADY!
        """
        super().__init__(visa_lock)

        # It is always good practice to initialize variables in the init

        # Instruments
        self.osa = None
        self.ps = None
        self.pm = None

        self.data = None

        if not self.check_necessary_instruments(instrument_list):
            raise ValueError("The necessary instruments for this experiment are not present!")

    def check_necessary_instruments(self, instrument_list):
        """
        Checks if the instruments necessary to perform the experiment are present.
        :param instrument_list: list of the available instruments
        :return: True if the necessary instruments are present, False otherwise.
        """

        # define parameters for Thorlabs pm.getSensorInfo()
        name = create_string_buffer(1024)  # sensor part number
        snr = create_string_buffer(1024)  # serial number
        message = create_string_buffer(1024)  # mystery date
        pType = c_int16()
        pStype = c_int16()
        pFlags = c_int16()

        for instr in instrument_list:
            if isinstance(instr, HP70951B):
                self.osa = instr
            if isinstance(instr, Keithley2635A):
                self.ps = instr
            if isinstance(instr, ThorlabsPowerMeter):
                instr.getSensorInfo(name, snr, message, byref(pType), byref(pStype), byref(pFlags))
                if (c_char_p(name.raw).value == SENSOR_1) or (c_char_p(name.raw).value == SENSOR_2):
                    self.pm = instr

        if (self.osa is not None) and (self.ps is not None) and (self.pm is not None):
            return True
        else:
            return False

    def get_description(self):
        """
        Returns a string with a brief summary of the experiment.
        """
        return "Sweeps voltage (tunes VCSEL) while recording measured voltage, osa peak and power, and SHG power into csv file"

    def get_name(self):
        """
        Returns a string with the experiment name
        """
        return "SHG_bandwidth"

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

        meas_volt_list = []
        peak_wl_list = []
        peak_dBm_list = []
        power_list = []

        [_, osa_settings_list] = self.osa.get_osa_parameters()

        # Sweep power supply voltage and get osa trace
        for ind, volt in enumerate(voltage_list):

            print('Setting power supply to %.4f V...' % volt)
            # Set the voltage
            self.ps.set_voltage(volt)

            # Read actual voltage
            meas_volt = self.ps.measure_voltage()
            print('Power Supply voltage set to %0.4f V' % meas_volt)
            meas_volt_list.append(meas_volt) #[V]

            # Wait [s]
            time.sleep(0.5)

            # Read osa highest peak and store amplitude, wavelength (VCSEL wavelength/power)
            [pk_wl, pk_amp] = self.osa.get_peak_info()
            peak_wl_list.append(pk_wl)
            peak_dBm_list.append(pk_amp)

            # Set Power meter to corresponding SHG wavelength and read power
            SHG_wavelength = float(pk_wl)/2 * 1e9  # convert to nm
            self.pm.setWavelength(c_double(SHG_wavelength))
            time.sleep(0.3)
            [SHG_power, _] = self.pm.get_powers()
            power_list.append(SHG_power)

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
                writer.writerow(meas_volt_list)  #[V]
                writer.writerow(peak_wl_list)  #[m]
                writer.writerow(peak_dBm_list) #[dBm]
                writer.writerow(power_list)  #[W]

                # Save the parameters in a .mat file
                time_tuple = time.localtime()
                params_filename = "params_%s_%d-%d-%d_%d-%d-%d.mat" % (
                    filename,
                    time_tuple[0],
                    time_tuple[1],
                    time_tuple[2],
                    time_tuple[3],
                    time_tuple[4],
                    time_tuple[5])

                io.savemat(params_filename, {'osa_settings_list': osa_settings_list,
                                             'voltage_list': voltage_list
                                             })

        self.data = [peak_wl_list, power_list]

        return [peak_wl_list, power_list]

    def required_params(self):
        """
        Returns a list with the keys that need to be specified in the params dictionary, in order for
        a measurement to be performed
        """
        return ["voltage_list"]

    def plot_data(self, canvas_handle, data=None):

        if data is None:
            if self.data is not None:
                data = self.data
            else:
                raise ValueError('plot_data was called before performing the experiment or providing data')

        wavelength = data[0]
        SHG_power = data[1]

        plot_graph(x_data=wavelength, y_data=SHG_power, canvas_handle=canvas_handle, xlabel='Wavelength (nm)',
                   ylabel='Power (uW)', title='SHG Power', legend=None)


import sys
from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg

class Window(QtGui.QMainWindow):
    def __init__(self):
        super(Window, self).__init__()
        self.setGeometry(100, 100, 1000, 500)
        self.setWindowTitle("SHG Power")

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
        self.xlabel = self.p_power.setLabel('bottom', text='Pump Wavelength', units='m')
        self.ylabel = self.p_power.setLabel('left', text='SHG Power', units='W')
        self.layout.addWidget(self.p_power, 0, 0)

        self.p_power_handle = self.p_power.plot(pen=(1, 3))

        self.show()

    def plot(self, x, y):
        self.p_power_handle.setData(x,y)

if __name__ == '__main__':
    # ------------------------------------------------------------
    # SAFETY LIMITS
    i_limit = 1e-9  # current limit [A]

    # OTHER PARAMETERS
    device = 'Dev1'
    pump_laser = 'CW980'
    pump_power = 20.5 #mW

    # EXPERIMENT PARAMETERS
    init_voltage = 53   # [V]
    end_voltage = 60 # [V]
    increment = 0.05  # Voltage increment
    voltage_list = np.arange(init_voltage, end_voltage+increment, increment) #end_voltage+1 or will stop at end_voltage-1
    # ------------------------------------------------------------

    # INSTRUMENTS
    # ps = KeysightE36106B(current_limit=i_limit)
    ps = Keithley2635A(current_compliance=i_limit, voltage_compliance=65) #A, V
    osa = HP70951B()
    pm = ThorlabsPowerMeter()

    # Initialize instruments
    ps.initialize()
    osa.initialize()
    pm.initialize()
    pm.setPowerAutoRange(c_int16(1))  #enable autorange


    file_name = "SHG_BW_%s_%d-%dV_%s_%3.2mWf" % (
        device, init_voltage, end_voltage, pump_laser, pump_power)  # Filename where to save csv data

    # SET UP THE EXPERIMENT
    instr_list = [osa, ps, pm]
    params = {"voltage_list": voltage_list}
    exp = SHG_bandwidth(instr_list)

    # RUN IT
    [pump_wl, SHG_power] = exp.perform_experiment(params, filename=file_name)

    # CLOSE INSTRUMENTS
    osa.close()
    ps.close()
    pm.close()

    # PLOT DATA
    app = QtGui.QApplication(sys.argv)
    GUI = Window()
    GUI.plot(pump_wl, SHG_power)
    sys.exit(app.exec_())


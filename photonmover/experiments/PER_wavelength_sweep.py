from photonmover.Interfaces.Experiment import Experiment
from photonmover.utils.plot_utils import plot_graph

# Interfaces/instruments necessary for the experiment
# - You use an Interface if any instrument of that category can be used
# - You use a specific instrument if you can only use that specific model
from photonmover.instruments.Power_meters.Thorlabs import ThorlabsPowerMeter
# from photonmover.instruments.Power_Supplies.KeysightE36106B import KeysightE36106B
from photonmover.instruments.Source_meters.Keithley2635A import Keithley2635A

# General imports
import time
import numpy as np
import csv

from ctypes import create_string_buffer, c_int16, c_double, byref, c_char_p, c_bool

# Power meter sensors
XPOL_SENSOR = b'S122C'  # 800-1700nm
PPOL_SENSOR = b'S122C'  # 800-1700nm


class PER_wavelength_sweep(Experiment):

    def __init__(self, instrument_list, visa_lock=None):
        """
        :param instrument_list: list of available instruments. IMPORTANT: WE ASSUME THAT THE INSTRUMENTS HAVE BEEN INITIALIZED ALREADY!
        """
        super().__init__(visa_lock)

        # It is always good practice to initialize variables in the init

        # Instruments
        self.pm_ppol = None
        self.pm_xpol = None
        self.ps = None

        self.data = None

        if not self.check_necessary_instruments(instrument_list):
            raise ValueError("The necessary instruments for this experiment are not present!")

    def check_necessary_instruments(self, instrument_list):
        """
        #pump_sensor=None, VCSEL_sensor=None
        Checks if the instruments necessary to perform the experiment are present.
        :param instrument_list: list of the available instruments, params - contains sensor info for assigning pump and VCSEL power meters
        :return: True if the necessary instruments are present, False otherwise.
        """

        # define parameters for Thorlabs pm.getSensorInfo()
        name = create_string_buffer(1024)  #sensor part number
        snr = create_string_buffer(1024)  #serial number
        message = create_string_buffer(1024)  #mystery date
        pType = c_int16()
        pStype = c_int16()
        pFlags = c_int16()

        # set flags to monitor if pump or VCSEL power meter assigned (for case where both power meters have same sensor)
        pm_xpol_set = 0
        pm_ppol_set = 0

        for instr in instrument_list:
            if isinstance(instr, ThorlabsPowerMeter):
                instr.getSensorInfo(name, snr, message, byref(pType), byref(pStype), byref(pFlags))
                if (c_char_p(name.raw).value == PPOL_SENSOR) and not pm_ppol_set:
                    self.pm_ppol = instr
                    pm_ppol_set = 1
                elif (c_char_p(name.raw).value == XPOL_SENSOR) and not pm_xpol_set:
                    self.pm_xpol = instr
                    pm_xpol_set = 1
            # if isinstance(instr, KeysightE36106B):
            if isinstance(instr, Keithley2635A):
                self.ps = instr

        if (self.pm_ppol is not None) and (self.pm_xpol is not None) and (self.ps is not None):
            return True
        else:
            return False

    def get_description(self):
        """
        Returns a string with a brief summary of the experiment.
        """
        return "Performs a PER measurement: sweeps VCSEL tuning voltage while measures VCSEL power thru polarizer.  Pump power is monitored with a 1% tap to a Thorlabs power meter."

    def get_name(self):
        """
        Returns a string with the experiment name
        """
        return "PER_wavelength_sweep"
        
    def perform_experiment(self, params, filename=None):
        """
        Performs the experiment, and saves the relevant data (if there is any)
        to the specified file (if given)
        :param params: dictionary of the parameters necessary for the experiment.
        :param filename: if specified, the data is saved in the specified file.
        :return: [pump_power_list, VCSEL_power_list]
        """

        params = self.check_all_params(params)

        VCSEL_wavelength = params["VCSEL_wavelength"]
        voltage = params["voltage"]

        # Set power meter wavelengths
        attribute = c_int16(0) #set wavelength
        meas_xpol_wavelength = c_double()
        meas_ppol_wavelength = c_double()

        self.pm_xpol.getWavelength(attribute, byref(meas_xpol_wavelength))
        self.pm_ppol.getWavelength(attribute, byref(meas_ppol_wavelength))
        if meas_ppol_wavelength is not VCSEL_wavelength:
            self.pm_ppol.setWavelength(c_double(VCSEL_wavelength))
        if meas_xpol_wavelength is not VCSEL_wavelength:
            self.pm_xpol.setWavelength(c_double(VCSEL_wavelength))

        print('Measured VCSEL pm wavelength is %d nm' % meas_xpol_wavelength.value)
        ppol_power_list = []
        xpol_power_list = []
        meas_volt_list = []

        # Sweep VOA voltage and get power
        for volt in voltage:

            print('Setting power supply voltage to %.4f V...' % volt)
            # Set the voltage
            self.ps.set_voltage(volt)

            meas_volt = self.ps.measure_voltage()
            print('Voltage set to %0.4f V' % meas_volt)
            meas_volt_list.append(meas_volt) #[V]

            # Wait s
            time.sleep(0.5)

            # Get VCSEL power at cross and parallel polarizations
            [ppol_power, _] = self.pm_ppol.get_powers()
            ppol_power_list.append(ppol_power) #[W]
            print('Measured PPol power is %.6f mW' % (ppol_power*1e3))

            [xpol_power, _] = self.pm_xpol.get_powers()
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
                writer.writerow(ppol_power_list)  #[W]
                writer.writerow(xpol_power_list)  #[W]
                writer.writerow(PER) #[dB]

        self.data = [meas_volt_list, PER]

        return [meas_volt_list, PER]
    
    def required_params(self):
        """
        Returns a list with the keys that need to be specified in the params dictionary, in order for
        a measurement to be performed
        """
        return ["xpol_sensor", "ppol_sensor", "VCSEL_wavelength", "voltage"]

    def plot_data(self, canvas_handle, data=None):
        
        if data is None:
            if self.data is not None:
                data = self.data
            else:
                raise ValueError('plot_data was called before performing the experiment or providing data')
        
        volt_list = data[0]
        PER = data[1]

        plot_graph(x_data=volt_list, y_data=PER, canvas_handle=canvas_handle, xlabel='Tuning Voltage (V)', ylabel='PER (dB)', title='PER sweep', legend=None)




import sys
from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg

class Window(QtGui.QMainWindow):
    def __init__(self):
        super(Window, self).__init__()
        self.setGeometry(100, 100, 1000, 500)
        self.setWindowTitle("PER Sweep")

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
        self.xlabel = self.p_power.setLabel('bottom', text='Tuning Voltage', units='V')
        self.ylabel = self.p_power.setLabel('left', text='PER', units='dB')
        self.layout.addWidget(self.p_power, 0, 0)

        self.p_power_handle = self.p_power.plot(pen=(1, 3))

        self.show()


    def plot(self, x, y):
        self.p_power_handle.setData(x,y)

if __name__ == '__main__':

    # ------------------------------------------------------------
    # SAFETY LIMITS
    i_limit = 0.003 #current limit

    # POWER METER SETTINGS
    # Check that power meter sensors (declared as global variables up top) are accurate
    # Note: pump sensor and VCSEL sensor can be the same type
    pump_wavelength = 976 #nm
    VCSEL_wavelength = 1252  #nm

    #OTHER DEVICE PARAMETERS
    device ='bf_PolStable'

    # EXPERIMENT PARAMETERS
    pump_power = 12.6 # [mW]
    IL = 0.78
    start_voltage = 0 # [V] Start tuning voltage in sweep
    end_voltage = -40  # [V] End tuning voltage in sweep
    increment = 1 # [V]
    if end_voltage > 0:
        volt_list = np.arange(start_voltage, end_voltage+1, increment)
    elif end_voltage <0:
        volt_list = np.arange(start_voltage, end_voltage-1, -increment)
    print(volt_list)
    # ------------------------------------------------------------

    # INSTRUMENTS
    # ps = KeysightE36106B(current_limit=i_limit)
    ps = Keithley2635A(current_compliance=0.002, voltage_compliance=151)
    pm1 = ThorlabsPowerMeter()
    pm2 = ThorlabsPowerMeter()

    # Set up power supply
    ps.initialize()

    # Set up power meters
    pm1.initialize(deviceIndex=0)
    pm2.initialize(deviceIndex=1)
    pm1.setPowerAutoRange(c_int16(1))  #enable autorange
    pm2.setPowerAutoRange(c_int16(1)) #enable autorange

    file_name = "PER_%s_%d-%dV_pm%d_OELand%dnm_%2.3fmW" % (device, start_voltage, end_voltage, VCSEL_wavelength, pump_wavelength, pump_power*IL) # Filename where to save csv data

    # SET UP THE EXPERIMENT
    instr_list = [pm1, pm2, ps]
    params = {"xpol_sensor": XPOL_SENSOR, "ppol_sensor": PPOL_SENSOR, "VCSEL_wavelength": VCSEL_wavelength, "voltage": volt_list}
    exp = PER_wavelength_sweep(instr_list)

    # RUN IT
    [volt_list, PER_list] = exp.perform_experiment(params, filename=file_name)


    # PLOT DATA
    app = QtGui.QApplication(sys.argv)
    GUI = Window()
    GUI.plot(volt_list, PER_list)
    sys.exit(app.exec_())

    # CLOSE INSTRUMENTS
    pm1.close()
    pm2.close()
    ps.close()
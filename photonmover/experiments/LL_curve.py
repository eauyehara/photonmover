from photonmover.Interfaces.Experiment import Experiment
from photonmover.utils.plot_utils import plot_graph

# Interfaces/instruments necessary for the experiment
# - You use an Interface if any instrument of that category can be used
# - You use a specific instrument if you can only use that specific model
from photonmover.instruments.Power_meters.Thorlabs import ThorlabsPowerMeter
from photonmover.instruments.Power_Supplies.AgilentE3633A import AgilentE3633A

# General imports
import time
import numpy as np
import csv

from ctypes import create_string_buffer, c_int16, c_double, byref, c_char_p, c_bool

# Power meter sensors
PUMP_SENSOR = b'S150C'  # 350-1100nm
VCSEL_SENSOR = b'S122C' #b'S155C'  # 800-1700nm


class LL_curve(Experiment):

    def __init__(self, instrument_list, visa_lock=None):
        """
        :param instrument_list: list of available instruments. IMPORTANT: WE ASSUME THAT THE INSTRUMENTS HAVE BEEN INITIALIZED ALREADY!
        """
        super().__init__(visa_lock)

        # It is always good practice to initialize variables in the init

        # Instruments
        self.pm_pump = None
        self.pm_VCSEL = None
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
        pm_pump_set = 0
        pm_VCSEL_set = 0

        for instr in instrument_list:
            if isinstance(instr, ThorlabsPowerMeter):
                instr.getSensorInfo(name, snr, message, byref(pType), byref(pStype), byref(pFlags))
                if (c_char_p(name.raw).value == PUMP_SENSOR) and not pm_pump_set:
                    self.pm_pump = instr
                    pm_pump_set = 1
                elif (c_char_p(name.raw).value == VCSEL_SENSOR) and not pm_VCSEL_set:
                    self.pm_VCSEL = instr
                    pm_VCSEL_set = 1
            if isinstance(instr, AgilentE3633A):
                self.ps = instr

        if (self.pm_pump is not None) and (self.pm_VCSEL is not None) and (self.ps is not None):
            return True
        else:
            return False

    def get_description(self):
        """
        Returns a string with a brief summary of the experiment.
        """
        return "Performs an LL measuremet: sweeps pump power by tuning VOA transmission with a power supply and measures VCSEL power.  Pump power is measured with a 1% tap to a Thorlabs power meter."

    def get_name(self):
        """
        Returns a string with the experiment name
        """
        return "LL_curve"
        
    def perform_experiment(self, params, filename=None):
        """
        Performs the experiment, and saves the relevant data (if there is any)
        to the specified file (if given)
        :param params: dictionary of the parameters necessary for the experiment.
        :param filename: if specified, the data is saved in the specified file.
        :return: [pump_power_list, VCSEL_power_list]
        """

        params = self.check_all_params(params)

        pump_wavelength = params["pump_wavelength"]
        VCSEL_wavelength = params["VCSEL_wavelength"]
        voltage = params["voltage"]
        IL = params["IL"]

        # Set power meter wavelengths
        attribute = c_int16(0) #set wavelength
        meas_pump_wavelength = c_double()
        meas_VCSEL_wavelength = c_double()

        self.pm_pump.getWavelength(attribute, byref(meas_pump_wavelength))
        self.pm_VCSEL.getWavelength(attribute, byref(meas_VCSEL_wavelength))
        if meas_pump_wavelength is not pump_wavelength:
            self.pm_pump.setWavelength(c_double(pump_wavelength))
        if meas_VCSEL_wavelength is not VCSEL_wavelength:
            self.pm_VCSEL.setWavelength(c_double(VCSEL_wavelength))

        print('Measured VCSEL pm wavelength is %d nm' % meas_VCSEL_wavelength.value)
        pump_power_list = []
        VCSEL_power_list = []
        meas_volt_list = []

        # Sweep VOA voltage and get power
        for volt in voltage:

            print('Setting VOA to %.4f V...' % volt)
            # Set the voltage
            self.ps.set_voltage(volt)

            meas_volt = self.ps.measure_voltage()
            print('VOA voltage set to %0.4f V' % meas_volt)
            meas_volt_list.append(meas_volt) #[V]

            meas_curr = self.ps.measure_current()
            print('VOA current is %0.4f A' % meas_curr)

            # Wait s
            time.sleep(0.5)

            #Read pump power through 1% tap, scale to power in 99% tap
            [pump_power_tap, _] = self.pm_pump.get_powers()
            pump_power = pump_power_tap*99*IL  #scale to full power
            pump_power_list.append(pump_power)  # [W]
            print('Pump power at %0.6f mW' % (pump_power*1e3))

            # Get VCSEL power
            [VCSEL_power, _] = self.pm_VCSEL.get_powers()
            VCSEL_power_list.append(VCSEL_power) #[W]
            print('Measured VCSEL power is %.6f mW' % (VCSEL_power*1e3))


        print('Finished LL curve')
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
                writer.writerow(pump_power_list)  #[W]
                writer.writerow(VCSEL_power_list)  #[W]

        self.data = [pump_power_list, VCSEL_power_list]

        return [pump_power_list, VCSEL_power_list]
    
    def required_params(self):
        """
        Returns a list with the keys that need to be specified in the params dictionary, in order for
        a measurement to be performed
        """
        return ["pump_sensor", "VCSEL_sensor", "pump_wavelength", "VCSEL_wavelength", "voltage", "IL"]

    def plot_data(self, canvas_handle, data=None):
        
        if data is None:
            if self.data is not None:
                data = self.data
            else:
                raise ValueError('plot_data was called before performing the experiment or providing data')
        
        pump_power = data[0]
        VCSEL_power = data[1]

        plot_graph(x_data=pump_power, y_data=VCSEL_power, canvas_handle=canvas_handle, xlabel='Pump Power (mW)', ylabel='VCSEL Power (mW)', title='LL curve', legend=None)




import sys
from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg

class Window(QtGui.QMainWindow):
    def __init__(self):
        super(Window, self).__init__()
        self.setGeometry(100, 100, 1000, 500)
        self.setWindowTitle("LL Curve")

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
        self.xlabel = self.p_power.setLabel('bottom', text='Pump Power', units='W')
        self.ylabel = self.p_power.setLabel('left', text='VCSEL Power', units='W')
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
    pump_wavelength = 1040 #nm
    VCSEL_wavelength = 1343  #nm

    IL = 0.750 #insertion loss measured as (WDM 980/1310 output) / (1% tap)*99% - scales 1% tap output to actual input to VCSEL fiber

    #OTHER DEVICE PARAMETERS
    device = 'Dev2a'
    tuning_voltage = 51 # [V]

    # EXPERIMENT PARAMETERS
    init_voltage = 3.5 #4.0  # [V] Minimum transmission on VOA (Note: when set to 5V, AgilentE3633A momentarily exceeds current limit when turning output on)
    end_voltage = 1.6 #2.24 #1.653  # [V] Maximum transmission on VOA
    num_points = 250  # Number of points between init and end current
    volt_list = np.linspace(init_voltage, end_voltage, num_points)
    # ------------------------------------------------------------

    # INSTRUMENTS
    ps = AgilentE3633A(current_limit=i_limit)
    pm1 = ThorlabsPowerMeter()
    pm2 = ThorlabsPowerMeter()

    # Set up power supply
    ps.initialize()


    # Set up power meters
    pm1.initialize(deviceIndex=0)
    pm2.initialize(deviceIndex=1)
    pm1.setPowerAutoRange(c_int16(1))  #enable autorange
    pm2.setPowerAutoRange(c_int16(1)) #enable autorange

    #file_name = 'LL_dev1_50V_pm1258_Solstis980nm'  # Filename where to save csv data
    file_name = "LL_%s_%dV_pm%d_pump%d" % (device, tuning_voltage, VCSEL_wavelength, pump_wavelength) # Filename where to save csv data

    # SET UP THE EXPERIMENT
    instr_list = [pm1, pm2, ps]
    params = {"pump_sensor": PUMP_SENSOR, "VCSEL_sensor": VCSEL_SENSOR, "pump_wavelength": pump_wavelength, "VCSEL_wavelength": VCSEL_wavelength, "voltage": volt_list, "IL": IL}
    exp = LL_curve(instr_list)

    # RUN IT
    [pump_power_list, VCSEL_power_list] = exp.perform_experiment(params, filename=file_name)


    # PLOT DATA
    app = QtGui.QApplication(sys.argv)
    GUI = Window()
    GUI.plot(pump_power_list, VCSEL_power_list)
    sys.exit(app.exec_())

    # CLOSE INSTRUMENTS
    pm1.close()
    pm2.close()
    ps.close()


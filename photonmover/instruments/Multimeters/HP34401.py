import sys
sys.path.insert(0, '../..')
import pyvisa as visa
import numpy as np
import time
from photonmover.Interfaces.Instrument import Instrument
GPIB_ADDR = "GPIB0::20::INSTR"  # VISA address

class HP34401(Instrument):
    """
    Code for controlling Agilent E3633A DC Power Supply
    """

    def __init__(self):
        super().__init__()

        # It is good practice to initialize variables in init
        self.gpib = None


    def initialize(self):
        """
        Initializes the instrument
        :return:
        """
        print('Opening connection to HP34401 Multimeter')

        rm = visa.ResourceManager()
        try:
            self.gpib = rm.open_resource(GPIB_ADDR, read_termination='\n', write_termination='\n', timeout=10000)
        except:
            raise ValueError('Cannot connect to HP34401 Multimeter')

        self.init_function()

    def init_function(self):
        """
        Nothing to do here
        """

    def close(self):
        print('Disconnecting HP34401 Multimeter')
        self.gpib.close()

    def error(self):
        print(self.gpib.query('SYST:ERR?'))

    def beep(self):
        """
        Sound the beeper!
        """
        self.gpib.write('SYST:BEEP')

    def measure_voltage(self, AC_DC, range='DEF', resolution='DEF'):
        """
        Measures the input voltage
        :param AC_DC: specify AC or DC
        :param range: 'MIN', 'MAX', 'DEF' (autorange)
        :param resolution in V
        :return: measured voltage
        """
        if AC_DC == "DC":
            meas_v = float(self.gpib.query('MEAS:VOLT:DC? %s,%s' % (range, resolution)))
        elif AC_DC == "AC":
            meas_v = float(self.gpib.query('MEAS:VOLT:AC? %s,%s' % (range, resolution)))
        else:
            print('Invalid input. Specify "AC" or "DC"')
        return meas_v

    def configure_dc_voltage(self):
        """
        Configure for fast single dc voltage reading
        :param sample_interval: time interval between voltage measurements [s]
        """

        # Configure for a DC voltage reading for the default range and resolution (5 digit slow)
        self.gpib.write('CONF:VOLT:DC DEF, DEF')
        self.gpib.write('VOLT:DC:NPLC 1')
        self.gpib.write('ZERO:AUTO OFF')

        # Set the number of readings per trigger to 1
        self.gpib.write('SAMP:COUNT 1')

        # Set internal trigger mode to trigger immediately
        self.gpib.write('TRIG:SOUR IMM')

    def measure_dc_voltage(self):
        """
        For fast DC voltage measurements
        :return:
        """
        volt = float(self.gpib.query('READ?'))
        return volt

    def measure_current(self, AC_DC, range='DEF', resolution='DEF'):
        """
        Measures the input current
        :param AC_DC: specify AC or DC
        :param range: 'MIN', 'MAX', 'DEF' (autorange)
        :param resolution in A (Default 5-1/2 digit resolution)
        :return: measured current
        """
        if AC_DC == "DC":
            meas_i = float(self.gpib.query('MEAS:CURR:DC? %s,%s' % (range, resolution)))
        elif AC_DC == "AC":
            meas_i = float(self.gpib.query('MEAS:CURR:AC? %s,%s' % (range, resolution)))
        else:
            print('Invalid input. Specify "AC" or "DC"')
        return meas_i

    def measure_resistance(self, range='DEF', resolution='DEF'):
        """
        Make a 2-wire ohms measurement
        :param range: 'MIN', 'MAX', 'DEF' (autorange)
        :param resolution: 'MIN', 'MAX', 'DEF'
        :return: resistance [ohms]
        """
        meas_r = float(self.gpib.query('MEAS:RES? %s,%s' % (range, resolution)))
        return meas_r

    def measure_numAcq_voltages(self, numAcq, sample_interval):
        """
        Note: Does not work properly! - actual sample interval from specified trigger delay longer than specified
        Better to use configure_dc_voltage() followed by measure_voltage() with time.sleep() for sampling at a specified interval

        Read multiple voltages at a specified sampling interval, triggered immediately
        :param numAcq: number of readings to acquire
        :param sample_interval: time interval between voltage measurements [s]
        :return: volt_array [V] array of measured voltages
        """

        # Configure for a DC voltage reading for the default range and resolution (5 digit slow)
        self.gpib.write('CONF:VOLT:DC DEF, DEF')
        self.gpib.write('VOLT:DC:NPLC 1')
        self.gpib.write('ZERO:AUTO OFF')

        # Set the number of readings per trigger (MIN = 1, MAX = 50,000)
        if numAcq >= 1 and numAcq <= 50000:
            self.gpib.write('SAMP:COUNT %s' % numAcq)
        else:
            print("Invalid input. Specify a number of acquisitions between 1 and 50,000")

        # Set trigger delay between each sample (MIN = 0.0015s, MAX = 3600s) - 1.5ms is default for default 5 digit resolution
        if sample_interval >= 0.0015 and sample_interval <= 3600:
            self.gpib.write('TRIG:DEL %s' % sample_interval)
        else:
            print("Invalid input. Specify a sample interval between 0 and 3600 [s]")

        # Set internal trigger mode to trigger immediately
        self.gpib.write('TRIG:SOUR IMM')

        # Set trigger state from "idle" to "wait-for-trigger" state
        volt_array = self.gpib.query_ascii_values('READ?', container=np.array)

        return volt_array

    def query_trig_delay(self):
        delay = self.gpib.query('TRIG:DEL?')
        return delay

    def query_integration_time(self):
        # self.gpib.write('VOLT:DC:NPLC 1')
        integration_time = self.gpib.query('VOLT:DC:NPLC?')
        return integration_time




if __name__ == '__main__':
    driver = HP34401()
    driver.initialize()
    print(driver.measure_voltage('DC'))
    # print(driver.measure_resistance())
    driver.beep()
    driver.close()

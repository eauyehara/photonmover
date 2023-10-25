import sys
sys.path.insert(0, '../..')
import pyvisa as visa
import time
from photonmover.Interfaces.Instrument import Instrument

GPIB_ADDR = "USB0::0x1313::0x804A::M00406896::INSTR"  # VISA adress


class ITC4001(Instrument):
    """
    Code for controlling HP 54750A oscilloscope
    """

    def __init__(self):
        super().__init__()

        # It is good practice to initialize variables in init
        self.gpib = None
        self.TEC_on = 0
        self.LD_on = 0

    def initialize(self):
        """
        Initializes the instrument
        :return:
        """
        print('Opening connection to ITC4001 LD/TEC driver')

        rm = visa.ResourceManager()
        try:
            self.gpib = rm.open_resource(GPIB_ADDR, read_termination='\n', write_termination='\n', timeout=10000)
        except:
            raise ValueError('Cannot connect to ITC4001 LD/TEC driver')

    def close(self):
        print('Disconnecting ITC4001 LD/TEC driver')
        self.gpib.close()

    def error(self):
        print(self.gpib.query('SYST:ERR?'))

    def turn_TEC_on(self):
        """
        Turns the source on
        :return:
        """
        self.gpib.write(":OUTP2 ON")
        self.TEC_on = 1

    def turn_TEC_off(self):
        """
        Turns the source off
        :return:
        """
        self.gpib.write(":OUTP2 OFF")
        self.TEC_on = 0

    def turn_LD_on(self):
        """
        Turns the source on
        :return:
        """
        self.gpib.write(":OUTP ON")
        self.LD_on = 1

    def turn_LD_off(self):
        """
        Turns the source off
        :return:
        """
        self.gpib.write(":OUTP OFF")
        self.LD_on = 0

    def measure_voltage(self):
        """
        Measures the laser voltage
        :return: measured voltage
        """
        meas_v = float(self.gpib.query('MEAS:VOLT?'))
        return meas_v

    def measure_temperature(self):
        """
        Measures the TEC temperature
        :return: measured temperature
        """
        meas_T = float(self.gpib.query('MEAS:TEMP?'))
        return meas_T

    def measure_current(self):
        """
        Measures the laser current
        :return: measured current
        """
        meas_i = float(self.gpib.query('MEAS:CURR?'))
        return meas_i

    def set_current_limit(self, current_limit):
        """
        Sets the laser current limit
        :param current_limit:
        :return:
        """
        self.gpib.write('SOUR:CURR:LIM %.6E' % current_limit)

    def set_current(self, current):
        """
        Sets the specified laser current and turns LD on
        :param current:
        :return:
        """
        self.gpib.write(":SOUR:CURR %.6E" % current)

        if not self.LD_on:
            self.turn_LD_on()

    def set_temperature(self, temperature):
        """
        Turns TEC on and sets the TEC temperature, then waits until stabilized within specified error
        :param temperature:
        :return:
        """
        T_error = 0.005  # [C]

        if not self.TEC_on:
            self.turn_TEC_on()

        self.gpib.write(":SOUR2:TEMP %.6E" % temperature)
        temp_start = self.measure_temperature()
        time.sleep(2.0)
        temp_later = self.measure_temperature()
        time.sleep(2.0)

        while abs(temp_later - temp_start) > T_error:
            temp_start = self.measure_temperature()
            time.sleep(2.0)
            temp_later = self.measure_temperature()
            time.sleep(2.0)
            print('Nope')
        print('Temperature setpoint reached')


if __name__ == '__main__':
    driver = ITC4001()
    driver.initialize()
    driver.set_current_limit(0.705)
    driver.set_current(0.01)
    print(driver.measure_current())
    driver.close()

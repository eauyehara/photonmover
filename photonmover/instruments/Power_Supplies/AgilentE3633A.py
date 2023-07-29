import sys
sys.path.insert(0, '../..')
import pyvisa as visa
import time
from Interfaces.Instrument import Instrument

GPIB_ADDR = "GPIB0::5::INSTR"  # VISA adress
DEFAULT_CURRENT_LIMIT = 1.0  # Default current limit in A


class AgilentE3633A(Instrument):
    """
    Code for controlling Agilent E3633A DC Power Supply
    """

    def __init__(self, current_limit=DEFAULT_CURRENT_LIMIT):
        super().__init__()

        # It is good practice to initialize variables in init
        self.gpib = None
        self.OUTPUT_on = 0
        self.cur_limit = current_limit

    def initialize(self):
        """
        Initializes the instrument
        :return:
        """
        print('Opening connection to Agilent E3633A DC Power Supply')

        rm = visa.ResourceManager()
        try:
            self.gpib = rm.open_resource(GPIB_ADDR, read_termination='\n', write_termination='\n', timeout=5000)
        except:
            raise ValueError('Cannot connect to Agilent E3633A DC Power Supply')

        self.init_function()

    def init_function(self):
        """
        Initializes a current limit for constant voltage mode
        """
        self.set_current_limit(self.cur_limit)

    def close(self):
        print('Disconnecting Agilent E3633A DC Power Supply')
        self.gpib.close()

    def error(self):
        print(self.gpib.query('SYST:ERR?'))


    def turn_output_on(self):
        """
        Turns the output on
        :return:
        """
        self.gpib.write(":OUTP ON")
        self.OUTPUT_on = 1

    def turn_output_off(self):
        """
        Turns the output off
        :return:
        """
        self.gpib.write(":OUTP OFF")
        self.OUTPUT_on = 0

    def set_volt_range(self, volt_range):
        """
        Toggles voltage range between 8V and 20V mode
        :return:
        """
        if volt_range == 'P8V':
            self.gpib.write("VOLT:RANG %s" %volt_range)
            print('Setting Voltage Range to 8V')
        elif volt_range == 'P20V':
            self.gpib.write("VOLT:RANG %s" % volt_range)
            print('Setting Voltage Range to 20V')
        else:
            print('Not a valid voltage range. Enter "P8V" or "P20V"')

    def measure_volt_range(self):
        """
        Queries voltage range
        :return: volt_range
        """
        volt_range = self.gpib.query('VOLT:RANG?')
        return volt_range


    def measure_voltage(self):
        """
        Measures the output voltage
        :return: measured voltage
        """
        meas_v = float(self.gpib.query('MEAS:VOLT?'))
        return meas_v


    def measure_current(self):
        """
        Measures the output current
        :return: measured current
        """
        meas_i = float(self.gpib.query('MEAS:CURR?'))
        return meas_i


    def set_voltage(self, voltage):
        """
        Checks if input voltage value is within range, then sets the output voltage and turns output on
        :param voltage:
        :return:
        """
        if voltage >= 0 and voltage <= 20:
            self.gpib.write(":SOUR:VOLT %.6E" % voltage)
        else:
            print('Not a valid voltage entry. Enter a value between 0 and 20V')

        if not self.OUTPUT_on:
            self.turn_output_on()


    def set_current_limit(self, current_limit):
        """
        Turns output off (if on) and set the current value
        :param current:
        :return:
        """
        if self.OUTPUT_on:
            self.turn_output_off()
        self.gpib.write(":SOUR:CURR %.6E" % current_limit)

    def set_volt_protection(self, volt_level):
        """
        Sets voltage protection level
        :param volt_level: voltage [V]
        :return:
        """
        volt_range = self.measure_volt_range()

        if volt_range == 'P8V':
            if volt_level >= 0 and volt_level <= 8.24:
                self.gpib.write("VOLT:PROT:LEV  %.6E" % volt_level)
            else:
                print('Specified voltage level is not valid for 8V voltage range (0-8.24V)')
        elif volt_range == 'P20V':
            if volt_level >= 0 and volt_level <= 20.6:
                self.gpib.write("VOLT:PROT:LEV  %.6E" % volt_level)
            else:
                print('Specified voltage level is not valid for 20V voltage range (0-20.6V)')


    def set_curr_protection(self, curr_level):
        """
        Sets current protection level
        :param level: current [A]
        :return:
        """
        volt_range = self.measure_volt_range()

        if volt_range == 'P8V':
            if curr_level >= 0 and curr_level <= 20:
                self.gpib.write("CURR:PROT:LEV  %.6E" % curr_level)
            else:
                print('Specified current level is not valid for 8V current range (0-20A)')
        elif volt_range == 'P20V':
            if curr_level >= 0 and curr_level <= 10:
                self.gpib.write("CURR:PROT:LEV  %.6E" % curr_level)
            else:
                print('Specified current level is not valid for 20V current range (0-10A)')

    def meas_volt_protection(self):
        """
        Queries voltage protection level
        :return:
        """
        volt_level = float(self.gpib.query('VOLT:PROT:LEV?'))
        return volt_level

    def meas_curr_protection(self):
        """
        Queries current protection level
        :return:
        """
        curr_level = float(self.gpib.query('CURR:PROT:LEV?'))
        return curr_level

    def clear_overvolt_condition(self):
        """
        Clears the overvoltage tripped circuit
        :return:
        """
        self.gpib.write("VOLT:PROT:CLE")

    def clear_overcurrent_condition(self):
        """
        Clears the overcurrent tripped circuit
        :return:
        """
        self.gpib.write("CURR:PROT:CLE")

if __name__ == '__main__':
    driver = AgilentE3633A(current_limit=0.003)
    driver.initialize()
    driver.set_current_limit(0.010)
    # print(driver.measure_current())

    driver.set_voltage(4.0)

    # print(driver.measure_current())
    #print(driver.measure_voltage())

    driver.close()

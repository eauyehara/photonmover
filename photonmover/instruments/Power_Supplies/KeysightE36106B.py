import sys
sys.path.insert(0, '../..')
import pyvisa as visa
import time
from Interfaces.Instrument import Instrument

GPIB_ADDR = "USB0::0x2A8D::0x1902::MY61001521::INSTR"  # VISA adress
DEFAULT_CURRENT_LIMIT = 0.003  # Default current limit in A


class KeysightE36106B(Instrument):
    """
    Code for controlling Keysight E36106B 100V DC Power Supply
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
        print('Opening connection to Keysight E36106B DC Power Supply')

        rm = visa.ResourceManager()
        try:
            self.gpib = rm.open_resource(GPIB_ADDR, read_termination='\n', write_termination='\n', timeout=5000)
        except:
            raise ValueError('Cannot connect to Keysight E36106B DC Power Supply')

        self.init_function()

    def init_function(self):
        """
        Initializes a current limit for constant voltage mode
        """
        self.set_current_limit(self.cur_limit)

    def close(self):
        print('Disconnecting Keysight E36106B DC Power Supply')
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


    def measure_voltage(self):
        """
        Measures the output voltage
        :return: measured voltage
        """
        meas_v = float(self.gpib.query('MEAS:VOLT:DC?'))
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
        if voltage >= 0 and voltage <= 100:
            self.gpib.write(":SOUR:VOLT %.6E" % voltage)
        else:
            print('Not a valid voltage entry. Enter a value between 0 and 100V')

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


        if volt_level >= 0 and volt_level <= 100:
            self.gpib.write("VOLT:PROT:LEV  %.6E" % volt_level)
        else:
            print('Specified voltage level is not valid (0-100V)')



    def set_curr_protection(self, curr_level):
        """
        Sets current protection level
        :param level: current [A]
        :return:
        """

        if curr_level >= 0 and curr_level <= 0.4:
            self.gpib.write("CURR:PROT:LEV  %.6E" % curr_level)
        else:
            print('Specified current level is not valid (0-0.4A)')

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
    driver = KeysightE36106B(current_limit=0.003)
    driver.initialize()
    # driver.set_current_limit(0.010)
    # print(driver.measure_current())

    driver.set_voltage(4.5)
    time.sleep(1.0)
    print(driver.measure_current())
    print(driver.measure_voltage())

    driver.close()

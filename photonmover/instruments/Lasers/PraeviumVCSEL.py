import pyvisa as visa
import time
import numpy as np
from scipy import io
from photonmover.Interfaces.Laser import TunableLaser
from photonmover.Interfaces.Instrument import Instrument
from photonmover.instruments.Source_meters.Keithley2635A import Keithley2635A


SWEEP_DWELL_TIME = 1.0  # Time to sleep at each wavelength when we do a
# tx curve by setting each wavelength at a time (in s)


class PraeviumVCSEL(Instrument, TunableLaser):
    """
    Tunes Praevium VCSEL wavelength with sourcemeter according to wavvolt '.mat' file
    """
    def __init__(self, smu, wavvolt_file):
        super().__init__()

        # Load wavelength over voltage data from mat file
        full_filename = wavvolt_file + ".mat"
        a = io.loadmat(full_filename)

        #Instruments
        self.smu = smu

        #Parameters
        self.wav = None
        self.volt = None
        self.power = None
        self.wav_range = a['peak_interp'][0] * 1e9 #[nm]
        self.volt_range = a['volt_interp'][0]  #[V]
        self.res = self.wav_range[1] - self.wav_range[0] #[nm] Wavelength resolution

        self.sweep_dwell_time = SWEEP_DWELL_TIME #[s]

        self.current_compliance = 100e-9 #[A]
        self.voltage_compliance = 81 #[V]


    def initialize(self):
        """
        Initializes the SMU, set current and voltage compliance
        """
        self.smu.initialize()
        self.smu.set_current_compliance(self.current_compliance)
        self.smu.set_voltage_compliance(self.voltage_compliance)
        self.smu.turn_on()


    def close(self):
        """
        Closes the instrument
        """
        self.smu.close()


    def turn_off(self):
        """
        Turn light off
        :return:
        """
        pass

    def turn_on(self):
        """
        Turn light on
        :return:
        """
        pass

    def set_power(self, power):
        """
        Set the power to the specified value (in mW)
        :return:
        """
        pass


    def set_wavelength(self, wavelength):
        """
        Set the wavelength to the specified value (in nm) from specified wavvolt file
        Wavelength-voltage tuning characteristics will change with CW/pulsed pumping, VCSEL temperature,
        or if device snaps down!
        :return: Interpolated wavelength from wavvolt file, measured smu voltage
        """

        #Check if wavelength within tuning range
        if (wavelength < np.min(self.wav_range) or (wavelength > np.max(self.wav_range)) or (wavelength > self.wav_range[0] and wavelength < self.wav_range[-1])):
            print("Wavelength out of range. Skipping wavelength")
            self.wav = None
        else:
            ind = np.argmin(np.abs(wavelength - self.wav_range))
            print('Setting wavelength to %.4f nm' % self.wav_range[ind])
            volt = self.volt_range[ind]
            self.smu.set_voltage(volt)
            time.sleep(1)
            self.volt = self.smu.measure_voltage()
            self.wav = self.wav_range[ind]
        return [self.wav, self.volt]


    def get_state(self):
        """
        Returns a list with the following elements:
        1. The current wavelength
        2. The current power
        """
        return [self.wav, self.power]

    def generate_valid_sweep(self, start_wl, stop_wl, step_wl):
        """
        Checks start and stop wavelengths and generates array of valid wavelengths (omits middle gap)
        """
        wl_arr = np.arange(start_wl, stop_wl+1, step_wl)
        valid_wl_arr = []

        if step_wl < self.res:
            print("Wavelength step below wavelength resolution")
        else:
            for wl in wl_arr:
                # Wavelength out of range
                if (wl < np.min(self.wav_range) or (wl > np.max(self.wav_range)) or (
                        wl > self.wav_range[0] and wl < self.wav_range[-1])):
                    pass
                # Wavelength in range
                else:
                    valid_wl_arr.append(wl)
        return valid_wl_arr

    def continuous_sweep(self, start_wl=1260, stop_wl=1340, step_wl=1, dwell_time=1, num_sweeps=2):
        """
        Sweeps wavelength num_sweeps times for alignment
        """
        wl_arr = self.generate_valid_sweep(start_wl, stop_wl, step_wl)
        for i in range(num_sweeps):
            for j, wl in enumerate(wl_arr):
                self.set_wavelength(wl)
                time.sleep(dwell_time)


if __name__ == '__main__':
    smu = Keithley2635A(current_compliance=100e-9, voltage_compliance=81)
    wavvolt_file = "wavvolt_Dev1a_25C_CW976_5.0mW"
    myLaser = PraeviumVCSEL(smu, wavvolt_file)
    myLaser.initialize()
    # valid_wl_arr = myLaser.generate_valid_sweep(1240, 1360,1)
    # print(valid_wl_arr)
    [wav, volt] = myLaser.set_wavelength(1270)
    # print(volt)
    # myLaser.continuous_sweep(dwell_time=1, num_sweeps=5)
    myLaser.close()


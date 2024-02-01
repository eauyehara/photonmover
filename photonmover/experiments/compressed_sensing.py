#Acquires dataset for compressed sensing
# - sweeps wavelength, acquires camera image at each wavelength, reads power from laser tap, output photodetector

from photonmover.Interfaces.Experiment import Experiment
from photonmover.utils.plot_utils import plot_graph


# Interfaces/instruments/experiments necessary for the experiment
# - You use an Interface if any instrument of that category can be used
# - You use a specific instrument if you can only use that specific model
from photonmover.Interfaces.Laser import Laser

from photonmover.instruments.DAQ.NI_DAQ import NiDAQ
from photonmover.instruments.Cameras.Xenics import Xenics
from photonmover.instruments.Lasers.SantecTSL210F import SantecTSL210F
from photonmover.instruments.Optical_spectrum_analyzers.HP70951B import HP70951B
from photonmover.instruments.Lasers.HPLightWave import HPLightWave

# General imports
import time
import numpy as np
import csv
import sys
import os
from pyqtgraph.Qt import QtGui
import pyqtgraph as pg

parent_dir = "C:\\Users\\elise\\Github\\photonmover\\photonmover\\experiments"
class compressed_sensing(Experiment):
    def __init__(self, instrument_list, visa_lock=None):
        """
        :param instrument_list: list of available instruments. IMPORTANT: WE ASSUME THAT THE INSTRUMENTS HAVE BEEN INITIALIZED ALREADY!
        """
        super().__init__(visa_lock)

        # It is always good practice to initialize variables in the init

        # Instruments
        self.daq = None
        self.laser = None
        self.camera = None
        self.osa = None

        # Save the last data obtained when the experiment was performed (for plotting purposes)
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
            if isinstance(instr, NiDAQ):
                self.daq = instr

            # if isinstance(instr, SantecTSL210F):
            if isinstance(instr, HPLightWave):
                self.laser = instr

            if isinstance(instr, Xenics):
                self.camera = instr

            if isinstance(instr, HP70951B):
                self.osa = instr

        if (
                self.daq is not None) and (
                self.laser is not None) and (
                self.camera is not None) and (
                self.osa is not None):
            return True
        else:
            return False


    def get_description(self):
        """
        Returns a string with a brief summary of the experiment.
        """
        return "Acquires a compressed sensing dataset"


    def get_name(self):
        """
        Returns a string with the experiment name
        """
        return "Compressed Sensing"

    def check_wavelength_sweep(self, RBW=0.1, wl_arr=None, start_wl=1250, stop_wl=1350, step_wl=0.1, filename=None):
        """
        Sweeps wavelength and reads set wavelength on OSA (Note: moves OSA window in 40nm intervals, or 0.1nm steps are lost)
        Returns set wavelength, measured wavelength on OSA
        """
        # Override start_wl, stop_wl, step_wl if wl_arr provided
        if wl_arr is not None:
            start_wl = wl_arr[0]
            stop_wl = wl_arr[-1]
            step_wl = 'custom'
            wavelength_arr = wl_arr
        else:
            wavelength_arr = np.arange(start_wl, stop_wl + step_wl, step_wl)

        # Set initial 0SA parameters
        osa_window = 20 # [nm]
        osa_start = str(start_wl - 10) + ' NM'
        osa_stop = str(start_wl + osa_window + 10) + ' NM'
        res_bw = str(RBW) + ' NM'
        self.osa.set_wl_axis(start_wl=osa_start, end_wl=osa_stop)
        self.osa.set_acq_bandwidth(res_bw=res_bw)

        # Set laser power, turn on laser
        self.laser.set_power(power=power)
        self.laser.set_wavelength(wavelength_arr[0])
        self.laser.turn_on()
        time.sleep(1.0)

        # Set laser wavelength, read OSA peak
        osa_wl = []
        osa_pk_amp = []
        wind_max = wavelength_arr[0] + osa_window

        for i, wl in enumerate(wavelength_arr):
            if wl >= wind_max:  # Shift the osa window
                osa_start = str(wind_max - 10) + ' NM'
                osa_stop = str(wind_max + osa_window + 10) + ' NM'
                self.osa.set_wl_axis(start_wl=osa_start, end_wl=osa_stop)
                wind_max = wind_max + osa_window
                time.sleep(2.0)

            self.laser.set_wavelength(wl)
            time.sleep(0.5)
            [wl, amp] = self.osa.get_peak_info()
            osa_wl.append(wl)
            osa_pk_amp.append(amp)
            time.sleep(0.5)

        self.laser.turn_off()

        if filename is not None:
            full_filename = "wlSweep_%s_%d-%d_step%s" % (filename, start_wl, stop_wl, str(step_wl)) + ".csv"
            with open(full_filename, 'w+') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(wavelength_arr)
                writer.writerow(osa_wl)
                writer.writerow(osa_pk_amp)

        print('Finished wavelength sweep')
        print('---------------------------------')

        return [wavelength_arr, osa_wl]

    # def perform_experiment(self, params, filename=None):
    #     """
    #     Performs the experiment, and saves the relevant data (if there is any) to the specified file (if given)
    #     :param params: dictionary of the parameters necessary for the experiment.
    #     :param filename: if specified, the data is saved in the specified file.
    #     :return:
    #     """
    #     # create path for data save
    #     os.mkdir(full_path)
    #
    #     params = self.check_all_params(params)
    #
    #     wavelengths = params["wavelengths"]
    #     power = params["power"]
    #     num_avg = params["num_avg"] # daq, tap acquisitions to average at each wavelength
    #     step = wavelengths[1] - wavelengths[0]
    #
    #     # Set laser power, turn on laser
    #     self.laser.set_wavelength(wavelengths[0])
    #     self.laser.set_power(power=power)
    #     self.laser.turn_on()
    #     time.sleep(2.0)
    #
    #     # Initialize data lists
    #     tap_meas = []
    #     PD_meas = []
    #
    #     # Measure wavelength[0] for photodiode data
    #     daq_read = self.daq.read_data(num_avg)
    #     tap_meas.append(np.mean(daq_read[0]))
    #     PD_meas.append(np.mean(daq_read[1]))
    #
    #
    #     # Loop through wavelengths
    #     print('Starting wavelength sweep')
    #     print('-----------------------------')
    #     for i, wl in enumerate(wavelengths[1:]):
    #         # Set laser wavelength
    #         self.laser.set_wavelength(wavelength=wl)
    #         time.sleep(0.3)
    #
    #         # Create camera filename for data save
    #         cam_filename = "%s_%4.1f.png" % (filename, (wl-step))
    #
    #         # Acquire camera image
    #         self.camera.get_frame(cam_filename)
    #
    #         daq_read = self.daq.read_data(num_avg)
    #         tap_meas.append(np.mean(daq_read[0]))
    #         PD_meas.append(np.mean(daq_read[1]))
    #
    #     cam_filename = "%s_%4.1f.png" % (filename, wavelengths[-1])
    #     self.camera.get_frame(cam_filename)
    #
    #     self.laser.turn_off()
    #
    #     # Save data in csv file
    #     csv_filename = "%s_sweep%4.1f-%4.1fnm_avg%d.csv" % (filename, wavelengths[0], wavelengths[-1], num_avg)
    #
    #     with open(csv_filename, 'w+') as csvfile:
    #         writer = csv.writer(csvfile)
    #         writer.writerow(wavelengths)
    #         writer.writerow(tap_meas)
    #         writer.writerow(PD_meas)
    #
    #     # Assign experiment data
    #     self.data = [wavelengths, PD_meas]
    #
    #     print('Finished wavelength sweep')
    #     print('-----------------------------')
    #
    #     return [wavelengths, PD_meas]

    def perform_experiment(self, params, filename=None):
        """
        Performs the experiment, and saves the relevant data (if there is any) to the specified file (if given)
        :param params: dictionary of the parameters necessary for the experiment.
        :param filename: if specified, the data is saved in the specified file.
        :return:
        """
        # create path for data save
        os.mkdir(full_path)

        params = self.check_all_params(params)

        wavelengths = params["wavelengths"]
        power = params["power"]
        num_avg = params["num_avg"]  # daq, tap acquisitions to average at each wavelength
        step = wavelengths[1] - wavelengths[0]

        # Set laser power, turn on laser
        self.laser.set_wavelength(wavelengths[0])
        self.laser.set_power(power=power)
        self.laser.turn_on()
        time.sleep(2.0)

        # Initialize data lists
        tap_meas = []
        PD_meas = []

        # Loop through wavelengths
        print('Starting wavelength sweep')
        print('-----------------------------')
        for i, wl in enumerate(wavelengths):
            # Set laser wavelength
            self.laser.set_wavelength(wavelength=wl)
            time.sleep(1.3)

            # Create camera filename for data save
            cam_filename = "%s_%4.1f.png" % (filename, wl)

            # Acquire camera image
            self.camera.get_frame(cam_filename)

            daq_read = self.daq.read_data(num_avg)
            tap_meas.append(np.mean(daq_read[0]))
            PD_meas.append(np.mean(daq_read[1]))

        self.laser.turn_off()

        # Save data in csv file
        csv_filename = "%s_sweep%4.1f-%4.1fnm_avg%d.csv" % (filename, wavelengths[0], wavelengths[-1], num_avg)

        with open(csv_filename, 'w+') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(wavelengths)
            writer.writerow(tap_meas)
            writer.writerow(PD_meas)

        # Assign experiment data
        self.data = [wavelengths, PD_meas]

        print('Finished wavelength sweep')
        print('-----------------------------')

        return [wavelengths, PD_meas]

    def required_params(self):
        """
        Returns a list with the keys that need to be specified in the params dictionary, in order for
        a measurement to be performed
        """
        return ["full_path", "wavelengths", "power", "num_avg"]


    def plot_data(self, canvas_handle, data=None):

        if data is None:
            if self.data is not None:
                data = self.data
            else:
                raise ValueError('plot_data was called before performing the experiment or providing data')

        x_axis = data[0]
        y_axis = data[1]

        plot_graph(
            x_data=x_axis,
            y_data=y_axis,
            canvas_handle=canvas_handle,
            xlabel='Wavelength (nm)',
            ylabel='PD Voltage (V)',
            title='Wavelength Sweep',
            legend=None)


class plot_window(QtGui.QMainWindow):
    def __init__(self):
        super(plot_window, self).__init__()
        self.setGeometry(100, 100, 1000, 500)
        self.setWindowTitle("Wavelength Sweep")

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
        self.xlabel = self.p_power.setLabel('bottom', text='Wavelength', units='nm')
        self.ylabel = self.p_power.setLabel('left', text='PD Voltage', units='V')
        self.layout.addWidget(self.p_power, 0, 0)

        self.p_power_handle = self.p_power.plot(pen=(1, 3))

        self.show()

    def plot(self, x, y):
        self.p_power_handle.setData(x, y)


if __name__ == '__main__':
    # ------------------------------------------------------------
    # EXPERIMENT PARAMETERS
    num_avg = 100
    sample = 'TiO2_4_Samp6ApCl3_recode3'

    # ------------------------------------------------------------
    # LASER PARAMETERS
    init_wavelength = 1484 #1260  # [nm] # Range on Santec TSL210F is 1260 to 1350 (noisy 1340-1350); Range on HP81640A is 1510-1640nm
    end_wavelength = 1584 #1340  # [nm]
    wl_step = 0.1  # [nm] wavelength step size
    wavelengths = np.arange(init_wavelength, end_wavelength+wl_step, wl_step)
    power = 0.4  # [mW]
    pm_tap_ch = 1  # power meter channel

    # ------------------------------------------------------------
    # DAQ PARAMETERS
    daq_dev = 'Dev1'
    daq_tap_ch = 'ai20'
    daq_PD_ch = 'ai16'
    daq_tap_min_val = 0
    daq_tap_max_val = 5
    daq_PD_min_val = 0
    daq_PD_max_val = 3
    sampling_freq = 100 #[Hz]

    daq_ch_list = [daq_dev + '/' + daq_tap_ch, daq_dev + '/' + daq_PD_ch]
    min_val_list = [daq_tap_min_val, daq_PD_min_val]
    max_val_list = [daq_tap_max_val, daq_PD_max_val]

    # ------------------------------------------------------------
    # Initialize instruments
    daq = NiDAQ()
    # laser = SantecTSL210F()
    laser = HPLightWave(tap_channel=pm_tap_ch, rec_channel=3)
    camera = Xenics()
    osa = HP70951B()

    daq.initialize()
    laser.initialize()
    camera.initialize()
    osa.initialize()

    # Configure instruments
    daq.configure_nsampl_acq(input_channels=daq_ch_list, num_points=num_avg, max_sampling_freq=sampling_freq, min_vals=min_val_list, max_vals=max_val_list)

    # Directory/Filename for data save
    time_stamp = time.localtime()
    data_dir = sample + '_%d_%d_%d' % (time_stamp[0], time_stamp[1], time_stamp[2])
    full_path = os.path.join(parent_dir, data_dir)
    filename = os.path.join(full_path, "cs_%s" % (sample))
    wlsweep_filename = 'HPLightWave_6' #'SantecTSL210F_23'

    # SET UP THE EXPERIMENT
    instr_list = [daq, laser, camera, osa]
    params = {"full_path": full_path, "wavelengths": wavelengths, "num_avg": num_avg, "power": power}
    exp = compressed_sensing(instr_list)

    # RUN IT
    [wave, PD_meas] = exp.perform_experiment(params, filename=filename)

    # wl_arr = temp
    # [set_wl, osa_wl] = exp.check_wavelength_sweep(RBW=0.1, start_wl=1260, stop_wl=1340, step_wl=0.1, filename=wlsweep_filename)
    # [set_wl, osa_wl] = exp.check_wavelength_sweep(RBW=0.1, start_wl=1460, stop_wl=1580, step_wl=0.1, filename=wlsweep_filename)

    # CLOSE INSTRUMENTS
    daq.close()
    camera.close()
    laser.close()
    osa.close()

    # PLOT DATA
    if exp.data is not None:
        x_axis = wave #set_wl
        y_axis = PD_meas #osa_wl
        app = QtGui.QApplication(sys.argv)
        GUI = plot_window()
        GUI.plot(x_axis, y_axis)
        sys.exit(app.exec_())
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

# temp = [1260,1260.1,1260.2,1260.5,1260.7,1260.9,1261.2,1261.4,1261.6,1261.9,1262.1,1262.3,1262.6,1262.8,1263.1,1263.3,1263.5,1263.7,1264,1264.2,1264.4,1264.7,1264.9,1265.1,1265.4,1265.6,1265.9,1266.1,1266.3,1266.5,1266.8,1267,1267.3,1267.5,1267.7,1268,1268.2,1268.4,1268.7,1268.9,1269.1,1269.4,1269.6,1269.9,1270.1,1270.3,1270.6,1270.8,1271,1271.3,1271.5,1271.8,1272.1,1272.2,1272.5,1272.7,1272.9,1273.3,1273.5,1273.7,1274,1274.2,1274.4,1274.7,1274.9,1275.2,1275.4,1275.6,1275.9,1276.1,1276.3,1276.6,1276.8,1277,1277.3,1277.5,1277.8,1278.1,1278.2,1278.5,1278.7,1279,1279.2,1279.4,1279.7,1279.9,1280.1,1280.4,1280.6,1280.8,1281.1,1281.2,1281.5,1281.8,1281.9,1282.2,1282.4,1282.7,1283,1283.1,1283.4,1283.7,1283.8,1284.1,1284.3,1284.5,1284.9,1285,1285.2,1285.6,1285.7,1286,1286.2,1286.4,1286.8,1286.9,1287.2,1287.4,1287.6,1288,1288.1,1288.4,1288.7,1288.8,1289.2,1289.3,1289.5,1289.9,1290,1290.3,1290.7,1290.8,1291.2,1291.5,1291.6,1291.9,1292.4,1292.5,1292.7,1293.1,1293.5,1293.9,1294.3,1294.7,1295.1,1295.5,1295.9,1296.3,1296.7,1296.8,1297.1,1297.5,1297.9,1298,1298.3,1298.6,1299,1299.4,1299.8,1300.1,1300.3,1300.6,1301,1301.3,1301.8,1302.1,1302.5,1302.9,1303.3,1303.7,1304.1,1304.4,1304.5,1304.8,1305,1305.3,1305.7,1306.1,1306.2,1306.5,1306.8,1306.9,1307.3,1307.7,1308.1,1308.5,1308.6,1308.8,1309.3,1309.7,1309.8,1310,1310.5,1310.9,1311,1311.3,1311.7,1312.1,1312.2,1312.5,1312.9,1313.3,1313.6,1313.7,1314.1,1314.5,1314.9,1315.3,1315.8,1316.2,1316.3,1316.5,1316.9,1317.3,1317.4,1317.7,1318.2,1318.5,1318.6,1318.9,1319,1319.3,1319.7,1320.2,1320.4,1320.6,1320.9,1321.3,1321.4,1321.7,1322.1,1322.6,1323,1323.4,1323.5,1323.8,1324.2,1324.6,1324.7,1325,1325.4,1325.8,1325.9,1326.2,1326.7,1327.1,1327.3,1327.5,1327.6,1327.9,1328.3,1328.5,1328.7,1328.8,1329.2,1329.3,1329.6,1330,1330.4,1330.6,1330.8,1330.9,1331.2,1331.6,1331.8,1332,1332.1,1332.5,1332.6,1332.9,1333,1333.3,1333.7,1333.8,1334.1,1334.5,1334.7,1335,1335.4,1335.6,1335.8,1335.9,1336.2,1336.3,1336.6,1337,1337.1,1337.4,1337.5,1337.9,1338.3,1338.6,1338.7,1339.1,1339.2,1339.5,1339.9,1340.4,1340.5,1340.8,1341,1341.2,1341.3,1341.6,1341.7,1342,1342.5,1342.6,1342.9,1343,1343.3,1343.7,1343.9,1344.1,1344.2,1344.6,1344.7,1345,1345.5,1345.6,1345.9,1346,1346.3,1346.4,1346.5,1346.8,1347.1,1347.6,1347.7,1348,1348.4,1348.5,1348.8,1348.9,1349.4,1349.7,1349.8,1350,1350.2,1350.4,1350.6]
# temp = [1260,1260.1,1260.5,1260.7,1260.9,1261.2,1261.4,1261.6,1261.9,1262.1,1262.3,1262.6,1262.8,1263.1,1263.3,1263.5,1263.7,1264,1264.2,1264.4,1264.7,1264.9,1265.1,1265.4,1265.6,1265.9,1266.1,1266.3,1266.5,1266.8,1267,1267.3,1267.5,1267.7,1268,1268.2,1268.4,1268.7,1268.9,1269.1,1269.4,1269.6,1269.9,1270.1,1270.3,1270.6,1270.8,1271,1271.3,1271.5,1271.8,1272.1,1272.2,1272.5,1272.7,1272.9,1273.3,1273.5,1273.7,1274,1274.2,1274.4,1274.7,1274.9,1275.2,1275.4,1275.6,1275.9,1276.1,1276.3,1276.6,1276.8,1277.3,1277.5,1277.8,1278.1,1278.2,1278.5,1278.7,1279,1279.2,1279.4,1279.7,1279.9,1280.1,1280.4,1280.6,1280.8,1281.1,1281.2,1281.5,1281.8,1281.9,1282.2,1282.4,1282.7,1283,1283.1,1283.4,1283.7,1283.8,1284.1,1284.3,1284.5,1284.9,1285,1285.6,1285.7,1286,1286.2,1286.4,1286.8,1286.9,1287.2,1287.4,1287.6,1288,1288.1,1288.4,1288.7,1288.8,1289.2,1289.3,1289.5,1289.9,1290,1290.3,1290.7,1290.8,1291.2,1292.7,1293.5,1293.9,1294.3,1294.7,1295.1,1295.5,1295.9,1296.3,1296.7,1296.8,1297.1,1297.5,1297.9,1298.3,1298.6,1299,1299.4,1299.8,1300.1,1300.6,1301,1301.8,1302.5,1302.9,1303.3,1303.7,1304.1,1304.5,1305,1305.3,1305.7,1306.1,1306.5,1306.9,1307.3,1307.7,1308.1,1308.5,1309.3,1309.7,1310.5,1310.9,1311.3,1311.7,1312.1,1312.5,1312.9,1313.3,1313.7,1314.1,1314.5,1315.3,1316.2,1316.5,1316.9,1317.3,1317.7,1318.2,1318.5,1318.9,1319.3,1319.7,1320.2,1320.6,1320.9,1321.3,1322.1,1322.6,1323,1323.4,1323.8,1324.2,1324.6,1325,1325.8,1325.9,1326.7,1327.1,1327.5,1327.9,1328.3,1329.2,1329.6,1330,1330.4,1330.8,1331.6,1332.1,1332.5,1332.6,1332.9,1333.3,1333.7,1334.1,1334.5,1335,1335.4,1335.8,1336.6,1337,1337.9,1338.3,1338.7,1339.1,1339.5,1339.9,1340.4,1340.5,1340.8,1341.2,1341.3,1341.6,1341.7,1342,1342.5,1342.6,1342.9,1343,1343.3,1344.6,1344.7,1345,1345.5,1346.4,1346.5,1346.8,1347.1,1347.6,1347.7,1348,1348.5,1348.8,1348.9,1350,1350.2,1350.4,1350.6]
temp = [1260,1260.5,1260.7,1260.9,1261.2,1261.4,1261.6,1261.9,1262.1,1262.3,1262.6,1262.8,1263.1,1263.3,1263.5,1263.7,1264,1264.2,1264.4,1264.7,1264.9,1265.1,1265.4,1265.6,1265.9,1266.1,1266.3,1266.5,1266.8,1267,1267.3,1267.5,1267.7,1268.2,1268.4,1268.7,1268.9,1269.1,1269.4,1269.6,1269.9,1270.1,1270.3,1270.6,1270.8,1271,1271.3,1271.5,1271.8,1272.2,1272.5,1272.7,1272.9,1273.3,1273.5,1273.7,1274,1274.2,1274.4,1274.7,1274.9,1275.2,1275.4,1275.6,1275.9,1276.1,1276.6,1276.8,1277.3,1277.5,1277.8,1278.1,1278.2,1278.5,1278.7,1279,1279.2,1279.4,1279.7,1279.9,1280.1,1280.4,1280.6,1280.8,1281.1,1281.5,1281.8,1281.9,1282.7,1283,1283.1,1283.4,1283.7,1283.8,1284.1,1284.5,1284.9,1285,1285.6,1285.7,1286,1286.2,1286.8,1286.9,1287.2,1287.4,1287.6,1288,1288.1,1288.4,1288.8,1289.2,1289.9,1290,1290.3,1290.8,1293.5,1293.9,1295.1,1295.9,1296.3,1296.7,1296.8,1297.1,1297.5,1297.9,1298.3,1299,1299.4,1299.8,1300.6,1301,1301.8,1302.5,1302.9,1303.3,1303.7,1304.1,1304.5,1305,1305.3,1305.7,1306.1,1306.5,1306.9,1307.3,1307.7,1308.1,1308.5,1309.3,1309.7,1310.5,1310.9,1311.3,1311.7,1312.1,1312.5,1313.3,1313.7,1314.5,1316.2,1318.2,1318.5,1318.9,1319.3,1320.2,1320.6,1322.6,1323,1323.4,1323.8,1324.2,1324.6,1326.7,1327.1,1327.5,1329.2,1329.6,1330,1330.4,1330.8,1332.1,1332.5,1333.3,1333.7,1334.5,1335,1335.4,1335.8,1336.6,1337,1338.3,1338.7,1339.1,1340.4,1340.8,1342,1342.6,1342.9,1343,1346.4,1346.5,1346.8,1347.7,1348.5]

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
        num_avg = params["num_avg"] # daq, tap acquisitions to average at each wavelength

        # Set laser power, turn on laser
        self.laser.set_wavelength(wavelengths[0])
        self.laser.set_power(power=power)
        self.laser.turn_on()
        time.sleep(2.0)

        # Initialize data lists
        # wave_meas = []
        tap_meas = []
        PD_meas = []

        # Loop through wavelengths
        print('Starting wavelength sweep')
        print('-----------------------------')
        for i, wl in enumerate(wavelengths):
            # Set laser wavelength
            self.laser.set_wavelength(wavelength=wl)
            time.sleep(0.3)

            # Reset lists containing average of num_sample measurements per wavelength
            # tap_to_ave = []
            # PD_to_ave = []

            # Create camera filename for data save
            cam_filename = "%s_%4.1f.png" % (filename, wl)

            # Acquire camera image
            self.camera.get_frame(cam_filename)

            # # Read actual laser wavelength
            # [_, wave, _] = self.laser.get_state()
            # wave_meas.append(wave)

            # Acquire num_avg tap and PD measurements (DAQ), store average
            # for iter in range(num_avg):
            #     # Reading power directly from laser over gpib
            #     [power_tap, _] = self.laser.get_powers()
            #     tap_to_ave.append(power_tap)
            #
            #     # Read daq channels
            #     daq_read = self.daq.task.read(number_of_samples_per_channel=num_avg)
            #     PD_to_ave.append(daq_read[0])
            #
            #     time.sleep(0.1)

            # Acquire num_avg and PD measurements from DAQ, store average
            # daq_read = self.daq.task.read(number_of_samples_per_channel=num_avg)
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
    sample = 'TiO2_4_HP_SampleApOpen_0.1nm'

    # ------------------------------------------------------------
    # LASER PARAMETERS
    init_wavelength = 1484 #1260  # [nm] # Range on Santec TSL210F is 1260 to 1350 (noisy 1340-1350); Range on HP81640A is 1510-1640nm
    end_wavelength = 1584 #1340  # [nm]
    wl_step = 0.1  # [nm] wavelength step size
    wavelengths = np.arange(init_wavelength, end_wavelength+wl_step, wl_step)
    # wavelengths = temp
    power = 0.4  # [mW]
    pm_tap_ch = 1  # power meter channel

    # ------------------------------------------------------------
    # DAQ PARAMETERS
    daq_dev = 'Dev1'
    daq_tap_ch = 'ai20'
    daq_PD_ch = 'ai16'
    daq_tap_min_val = 0
    daq_tap_max_val = 1
    daq_PD_min_val = -1
    daq_PD_max_val = 1
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
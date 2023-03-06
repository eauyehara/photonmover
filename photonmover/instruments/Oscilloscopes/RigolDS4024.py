import sys
sys.path.insert(0, '../..')
import pyvisa as visa
import numpy as np
import time
import csv
from Interfaces.Instrument import Instrument


GPIB_ADDR = "USB0::0x1AB1::0x04B1::DS4A201100188::INSTR"  # VISA adress


class RigolDS4042(Instrument):
    """
    Code for controlling Rigol DS4042 Oscilloscope
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
        print('Opening connnection to Rigol Oscilloscope')

        rm = visa.ResourceManager()
        try:
            self.gpib = rm.open_resource(GPIB_ADDR, timeout=10000)
        except:
            raise ValueError('Cannot connect to Rigol Oscilloscope')


    def close(self):
        print('Disconnecting Rigol Oscilloscope')
        self.gpib.close()

    def autoscale(self):
        self.gpib.write(":AUT")

    def clear(self):
        """
        Clears all the waveforms in the display
        :return:
        """
        self.gpib.write(":CLE")

    def run(self):
        """
        Start the oscilloscope
        :return:
        """
        self.gpib.write(":RUN")

    def stop(self):
        """
        Stop the oscilloscope
        :return:
        """
        self.gpib.write(":STOP")

    def single(self):
        """
        Single acquisition
        :return:
        """
        self.gpib.write(":SING")

    def force_trigger(self):
        """
        Force a trigger
        :return:
        """
        self.gpib.write(":TFOR")

    def set_acq_averages(self, num_avg):
        """
        Sets averaging to the specified number
        :param num_avg:
        :return:
        """

        if num_avg > 1024:
            print("Specified averages are too high. Setting to maximum.")
            num_avg = 1024

        if num_avg % 2 != 0:
            print("The number of averages has to be a power of 2. Setting to the closest.")
            num_avg = round(np.log2(num_avg))

        self.gpib.write("ACQ:AVER %d" % num_avg)

    def set_acq_type(self, mode):
        """
        Sets the acquisition type for the oscilloscope
        :param mode: "NORM", "AVER", "PEAK", "HRES"
        :return:
        """

        if mode not in ["NORM", "AVER", "PEAK", "HRES"]:
            print("Acquisition mode not supported. Doing nothing.")
            return

        self.gpib.write(":ACQ:TYPE %s" % mode)

    def set_bw(self, channel, bw=0):
        """
        Sets the bandwidth limit for the specified channel
        :param channel: 1, 2, 3 or 4
        :param bw: 0 if OFF (default), 1 if 20MHz bw limit, 2 if 100MHz BW limit
        :return:
        """

        if channel not in [1, 2, 3 ,4]:
            print("Channel not correct. Doing nothing.")
            return

        if bw not in [0, 1, 2]:
            print("Bandwidth option not correct. Doing nothing")
            return

        if bw == 0:
            self.gpib.write(":CHAN%d:BWL OFF" % channel)

        if bw == 1:
            self.gpib.write(":CHAN%d:BWL 20M" % channel)

        if bw == 2:
            self.gpib.write(":CHAN%d:BWL 100M" % channel)

    def set_coupling_type(self, channel, coupling):
        """
        Set the coupling type (AC, DC or GND)
        :param channel: 1, 2, 3 or 4
        :param coupling: "AC", "DC", "GND"
        :return:
        """

        if channel not in [1, 2, 3 ,4]:
            print("Channel not correct. Doing nothing.")
            return

        if coupling not in ["AC", "DC", "GND"]:
            print("Coupling option not correct. Doing nothing")
            return

        self.gpib.write(":CHAN%d:COUP %s" % (channel, coupling))

    def channel_display(self, channel, on):
        """
        Turns on or off the display of the specified channel
        :param channel: 1, 2, 3 or 4
        :param on: 0 for off, 1 for on
        :return:
        """

        if channel not in [1, 2, 3 ,4]:
            print("Channel not correct. Doing nothing.")
            return

        if on not in [0, 1]:
            print("Display option not correct. Doing nothing")
            return

        self.gpib.write(":CHAN%d:DISP %d" % (channel, on))

    def set_impedance(self, channel, imped):
        """
        Sets the channel impedance
        :param channel: 1, 2, 3, or 4
        :param imped: OMEG for 1MOhm, FiFT for 50ohms
        :return:
        """

        if channel not in [1, 2, 3, 4]:
            print("Channel not correct. Doing nothing.")
            return

        if imped not in ["OMEG", "FIFT"]:
            print("Impedance option not correct. Doing nothing")
            return

        self.gpib.write(":CHAN%d:IMP %s" %(channel, imped))

    def set_vertical_range(self, channel, scale, offset):
        """
        Set the vertical range parameters for the specified channel. If they are
        None, nothing happens.
        :param channel: 1, 2, 3 or 4
        :param scale: for 50ohm impedance, 1x probe ratio: 1mV/div to 1 V/div; for 1Mohm impedance and 1x probe ratio: 1mV/div to 5V/div
        :param offset: see programming manual
        :return:
        """

        if channel not in [1, 2, 3 ,4]:
            print("Channel not correct. Doing nothing.")
            return

        if scale is not None:
            self.gpib.write(":CHAN%d:SCAL %.3f" % (channel, scale))

        if offset is not None:
            self.gpib.write(":CHAN%d:OFFS %.3f" % (channel, offset))

    def set_horizonal_scale(self, scale):
        """
        Sets the horizontal scale in s/div
        :param scale: 2ns to 1000s in Y-T mode, 200ms to 1000s in Roll mode
        :return:
        """
        self.gpib.write(":TIM:SCAL %.7f" % scale)

    def measure_item(self, channel, item):
        """
        Measures the specified item in the specified channel
        :param channel: 1, 2, 3 or 4
        :param item: "VMAX", "VMIN", "VPP", "VAVG", "PER", "FREQ"
        :return: the specified item value
        """

        if channel not in [1, 2, 3, 4]:
            print("Channel not correct. Doing nothing.")
            return

        if item not in ["VMAX", "VMIN", "VPP", "VAVG", "PER", "FREQ"]:
            print("Specified item to measure not correct. Doing nothing")
            return

        return float(self.gpib.query_ascii_values(":MEAS:%s CHAN%d" % (item, channel)))
        # time.sleep(1)
        # return float(self.gpib.query_ascii_values(":MEAS:ITEM? %s,CHAN%d" % (item, channel))[0])

    def set_trigger(self, mode, coupling, trig_sweep, channel, level):
        """
        Sets the trigger with the specified parameters.
        :param mode: Trigger mode. One of "EDGE", "PULSE", "RUNT", "WIND", "SLOPE", "NEDG", "PATT", "DEL"
        :param coupling: Coupling type. One of "AC", "DC", "LFR", "HFR"
        :param trig_number: Sets when the trigger fires. One of "AUTO", "NORM", "SING"
        :param channel: source of the trigger
        :param level: level for the trigger to fire
        :return:
        """

        if mode in ["EDGE", "PULSE", "RUNT", "NEDG", "SLOP", "VID","PATT", "RS232", "IIC", "SPI", "CAN", "FLEX", "USB"]:
            self.gpib.write(":TRIG:MODE %s" % mode)

        if coupling in ["AC", "DC", "LFR", "HFR"]:
            self.gpib.write(":TRIG:COUP %s" % coupling)

        if trig_sweep in ["AUTO", "NORM", "SING"]:
            self.gpib.write(":TRIG:SWE %s" % trig_sweep)

        if channel in [1, 2, 3, 4]:
            self.gpib.write(":TRIG:%s:SOUR CHAN%d" % (mode, channel))

        if level is not None:
            self.gpib.write(":TRIG:%s:LEV %.4f" % (mode, level))

    def read_waveform(self, channels, file_name=None):
        """
        Reads the waveform in the specified channels
        :param channels: list with the channels whose waveform we w ant to obtain.
        :param file_name: if specified, it will save the data with the specified file name. Do not include the ".csv". One file for each channel
        :return: 2 lists, each with n elemnts, where n is the number of channels.
                List 1: [preamble_channel1, preamble_channel2, ...]
                List 2: [channel1_data, channel2_data, ...]
        """

        all_preambles = []
        all_waveforms = []

        # Set waveform reading mode to normal - reads waveform displayed on screen
        # self.gpib.write(":STOP")
        self.gpib.write(":WAV:MODE NORM")

        # Set to send ascii
        self.gpib.write(":WAV:FORM ASCII")

        for c in channels:
            if c not in [1, 2, 3, 4]:
                print("Specified channel not correct. Skipping it")
                continue

            # Choose source
            self.gpib.write(":WAV:SOUR CHAN%d" % c)

            self.gpib.write("WAV:DATA?")
            data = self.gpib.read_raw()
            #print(data)
            raw_data = data[11:]
            #print(raw_data)

            wav_data_str = (str(raw_data)[2:-4]).split(',')

            wav_data = [float(i) for i in wav_data_str]
            #print(wav_data)

            preamble = self.gpib.query_ascii_values("WAV:PRE?")

            # Save the data if necessary. Each channel will be stored in a different file
            if file_name is not None:

                # Create the csv file
                file_name_chan = file_name + "_channel_" + str(c) + ".csv"

                with open(file_name_chan, 'w+') as csvfile:
                    writer = csv.writer(csvfile)
                    writer.writerow(preamble)
                    writer.writerow(wav_data)

            all_preambles.append(preamble)
            all_waveforms.append(wav_data)

            # self.gpib.write(":RUN")

        return all_preambles, all_waveforms

    def get_channel_states(self):
        """
        Checks which channels are activated.
        Returns list of booleans.
        """
        return [bool(self.gpib.query_ascii_values(":CHANnel{:d}:DISPlay?".format(x))[0]) for x in [1,2,3,4]]


    def set_memory_depth(self, memory_depth):
        """
        Sets memory depth
        - if 1 channel enabled: 14kPts, 140kPts, 1.4MPts, 14MPts, 140MPts
        - if 2 channels ([Ch1 and Ch2] or [Ch3 and Ch4]): 7kPts, 70kPts, 700kPts, 7MPts, 70MPts
        :param memory_depth:
        :return:
        """
        chans_enabled = self.get_channel_states()
        num_chans_enabled = np.sum(chans_enabled)

        self.gpib.write(":RUN")  # Can only set memory depth in Run state

        if num_chans_enabled == 1:
            assert memory_depth in ['AUTO', 14e3, 140e3, 1.4e6, 14e6, 140e6]
        else: # num_chans_enabled > 1 and num_chans_enabled <=4:
            all_channels = np.array(range(1, 5))
            if np.all(all_channels[chans_enabled] == np.array([1,2])) or np.all(all_channels[chans_enabled] == np.array([3,4])):
                assert memory_depth in ['AUTO', 7e3, 70e3, 700e3, 7e6, 70e6]
            else:
                assert memory_depth in ['AUTO', 14e3, 140e3, 1.4e6, 14e6, 140e6]
                #still acts as single channel

        if memory_depth == 'AUTO':
            self.gpib.write(":ACQ:MDEP %s" %memory_depth)
        else:
            self.gpib.write(":ACQ:MDEP %d" %memory_depth)


    def get_memory_depth(self):
        """
        Get memory depth [pts]
        :return: memory_depth
        """
        return int(self.gpib.query(":ACQ:MDEP?"))


    def get_status(self):
        """
        Get status of data read
        :return:
        """
        ans = str(self.gpib.query(":WAV:STAT?")).split(',')
        status = ans[0]
        pts_read = ans[1].strip('\n')
        return [status, pts_read]


    def read_FFT(self, file_name=None):
        """
        Reads the FFT spectrum on screen and corresponding channel waveform (assumes channel, span, center already optimized on front panel)
        :param file_name: if specified, it will save the data and params_FFT with the specified file name. Do not include the ".csv".
        :return: 2 lists, one with FFT parameters, one with the FFT data
                List 1: [params_FFT]
                List 2: [FFT_data]
        """
        # Set horizontal position of FFT to 0 Hz
        self.gpib.write("CALC:FFT:HOFF 0")

        # Set waveform reading mode to normal - reads waveform displayed on screen
        self.gpib.write(":WAV:MODE NORM")

        # Set to send ascii
        self.gpib.write(":WAV:FORM ASCII")

        # Choose FFT as source
        self.gpib.write(":WAV:SOUR FFT")

        # Read FFT data
        self.gpib.write("WAV:DATA?")
        data = self.gpib.read_raw()
        # print(data)
        raw_data = data[11:]
        # print(raw_data)

        #Parse raw data and convert to float
        FFT_data_str = (str(raw_data)[2:-4]).split(',')
        # print(FFT_data_str)
        FFT_data = [float(i) for i in FFT_data_str]

        # # Query FFT params
        HCEN = float(self.gpib.query("CALC:FFT:HCEN?")) #Center frequency [Hz])
        HOFFS = float(self.gpib.query("CALC:FFT:HOFF?")) #Horizontal position of FFT [Hz]
        HSC = float(self.gpib.query("CALC:FFT:HSC?") )#Horizontal coefficient of FFT operation [1, 2, 3, 4]
        HSPAN = float(self.gpib.query("CALC:FFT:HSP?")) #Horizontal scale (freq value per grid) [Hz/div]
        FFT_source = self.gpib.query("CALC:FFT:SOUR?") #source channel of FFT
        # FFT_window = self.gpib.query("CALC:FFT:WIND?") #{RECT,HANN,HAMM,BLACK}
        
        # VSMODE = self.gpib.query("CALC:FFT:VSM?") #vertical scale type [DB or VRMS]
        VSCALE = float(self.gpib.query("CALC:FFT:VSC?"))  # vertical scale of FFT [dBV/div (dBm/div in 50ohm) or V/div]
        VOFFSET = float(self.gpib.query("CALC:FFT:VOFF?"))  # vertical position of FFT [dBV (dBm in 50ohm) or V]
        FFT_params = [HCEN, HOFFS, HSC, HSPAN, VSCALE, VOFFSET]


        # Read waveform data of FFT source channel
        self.gpib.write(":WAV:SOUR %s" % FFT_source)
        self.gpib.write("WAV:DATA?")
        data = self.gpib.read_raw()
        raw_data = data[11:]

        wav_data_str = (str(raw_data)[2:-4]).split(',')
        wav_data = [float(i) for i in wav_data_str]

        preamble = self.gpib.query_ascii_values("WAV:PRE?")

        # Save the data if necessary. Each channel will be stored in a different file
        if file_name is not None:
            # Create the csv file
            file_name_chan = file_name + "_FFT.csv"

            with open(file_name_chan, 'w+') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(preamble)
                writer.writerow(wav_data)
                writer.writerow(FFT_params)
                writer.writerow(FFT_data)

        return FFT_params, FFT_data


    def read_internal_memory(self, channels, memory_depth=14e3, file_name=None):
        """
        Starts a new acquisition and reads the waveform data in internal memory
        :param channels: list with the channels whose waveform we w ant to obtain.
        :param file_name: if specified, it will save the data with the specified file name. Do not include the ".csv". One file for each channel
        :return: 2 lists, each with n elements, where n is the number of channels.
                List 1: [preamble_channel1, preamble_channel2, ...]
                List 2: [channel1_data, channel2_data, ...]
        """

        time_list = []
        all_waveforms = []
        wav_data = []

        # Set memory depth - raw mode grabs points 1 to current max memory depth
        self.set_memory_depth(memory_depth)

        #Scope needs to be in run state after last memory buffer read or buffer will be empty
        self.gpib.write(":RUN") #In case scope is still in STOP state from last read
        time.sleep(5) #if buffer is still empty, need to increase this time

        self.gpib.write(":STOP")  # Can only read internal memory when oscilloscope in Stop state
        self.gpib.write(":WAV:MODE RAW")  # Set waveform reading mode to raw
        self.gpib.write(":WAV:FORM ASCII")  # Set return format of waveform data to binary
        self.gpib.write(":WAV:POIN %d" % memory_depth) # set number of points read to current memory depth

        for c in channels:
            if c not in [1, 2, 3, 4]:
                print("Specified channel not correct. Skipping it")
                continue

            #Choose source
            self.gpib.write(":WAV:SOUR CHAN%d" % c)

            #Get time axis data
            t_orig = float(self.gpib.query_ascii_values(":WAV:XOR?")[0])
            t_inc = float(self.gpib.query_ascii_values(":WAV:XINC?")[0])
            t_ref = float(self.gpib.query_ascii_values(":WAV:XREF?")[0])
            time_data = np.arange(0, memory_depth*t_inc, t_inc) + t_orig + t_ref

            self.gpib.write(":WAV:RES")
            self.gpib.write(":WAV:BEG")

            # Collect data until status is set to IDLE
            [status, _] = self.get_status()
            while status == "READ":
                print(status)
                self.gpib.write(":WAV:DATA?")
                data = self.gpib.read_raw()
                if data == b'#9000000000\n':
                    continue

                parse_data = str(data[11:])[2:-4].split(',')
                wav_data += [float(i) for i in parse_data]
                time.sleep(0.1)
                [status, _] = self.get_status()

            if status == "IDLE":
                print(status)
                self.gpib.write(":WAV:DATA?")
                data = self.gpib.read_raw()

                if data == b'#9000000000\n':
                    print('Buffer empty')
                else:
                    parse_data = str(data[11:])[2:-4].split(',')
                    wav_data += [float(i) for i in parse_data]
                self.gpib.write(":WAV:END")

            print("%d points read" % len(wav_data))

            # preamble = self.gpib.query_ascii_values("WAV:PRE?")

            # Save the data if necessary. Each channel will be stored in a different file
            if file_name is not None:
                # Create the csv file
                file_name_chan = file_name + "_channel_" + str(c) + ".csv"

                with open(file_name_chan, 'w+') as csvfile:
                    writer = csv.writer(csvfile)
                    writer.writerow(time_data)
                    writer.writerow(wav_data)

            time_list.append(time_data)
            all_waveforms.append(wav_data)

        self.gpib.write(":RUN")

        return time_list, all_waveforms



if __name__ == '__main__':

    osc = RigolDS4042()
    osc.initialize()

    # osc.read_FFT('Rigol_ch1_PDB210A2_sigMon_976CW_span874kHz_RECT_10mVdiv_ACcoupled')
    # osc.read_waveform([1], 'Rigol_PolStable_0V_ET3500F_polController_GTPol_parallel_OEland5mW')
    osc.read_internal_memory([1], memory_depth=14e3, file_name='Rigol_PolStable_0V_150MHzInGaAs_orthog_OEland5mW_memdep14e3_20nsdiv')

    osc.close()

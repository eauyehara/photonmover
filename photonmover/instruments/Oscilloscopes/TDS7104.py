import pyvisa as visa
import numpy as np
import time
import csv
from photonmover.Interfaces.Instrument import Instrument
import binascii

GPIB_ADDR = "GPIB0::15::INSTR"  # VISA adress


class TDS7104(Instrument):
    """
    Code for controlling Tektronix TDS7104 Oscilloscope
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
        print('Opening connnection to Tektronix Oscilloscope')

        rm = visa.ResourceManager()
        try:
            self.gpib = rm.open_resource(GPIB_ADDR, timeout=20000)
        except:
            raise ValueError('Cannot connect to Tektronix Oscilloscope')


    def close(self):
        print('Disconnecting Tektronix Oscilloscope')
        self.gpib.close()

    def autoset(self):
        self.gpib.write("AUTOS")

    def run(self):
        """
        Start the oscilloscope
        :return:
        """
        self.gpib.write("ACQ:STATE:RUN")

    def stop(self):
        """
        Stop the oscilloscope
        :return:
        """
        self.gpib.write("ACQ:STATE:STOP")

    def single(self):
        """
        Single acquisition
        :return:
        """
        self.gpib.write("ACQ:STOPA SEQ")

    def force_trigger(self):
        """
        Force a trigger
        :return:
        """
        self.gpib.write("TRIG FORC")

    def set_acq_averages(self, num_avg):
        """
        Sets Acquisition mode to average and averages for the specified number
        :param num_avg:
        :return:
        """
        self.set_acq_mode("AVE")

        if num_avg > 10000:
            print("Specified number of averages too high. Setting to maximum 10,000.")
            num_avg = 10000

        if num_avg < 2:
            print("Specified number of averages too low. Setting to minimum 2.")
            num_avg = 2

        self.gpib.write("ACQ:NUMAV %d" % num_avg)

    def set_acq_mode(self, mode):
        """
        Sets the acquisition type for the oscilloscope
        :param mode: {SAMple|PEAKdetect|HIRes|AVErage|ENVelope|WFMDB}
        :return:
        """

        if mode not in ["SAM", "SAMple","PEAK", "PEAKdetect", "HIR", "HIRes","AVE","AVErage","ENV","ENVelope","WFMDB"]:
            print("Acquisition mode not supported. Doing nothing.")
            return

        self.gpib.write("ACQ:MOD %s" % mode)

    def set_bw(self, channel, bw="FUL"):
        """
        Sets the bandwidth limit for the specified channel
        :param channel: 1, 2, 3 or 4
        :param bw: ["TWEnty", "TWOfifty", "FULl"]
        :return:
        """

        if channel not in [1, 2, 3 ,4]:
            print("Channel not correct. Doing nothing.")
            return

        if bw not in ["TWE", "TWEnty", "TWO", "TWOfifty", "FUL", "FULl"]:
            print("Bandwidth option not correct. Doing nothing")
            return
        self.gpib.write("CH%d:BAN %s" %(channel, bw))


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

        self.gpib.write("CH%d:COUP %s" % (channel, coupling))

    def channel_on(self, channel, on):
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
        if on == 1:
            set = "ON"
        elif on == 0:
            set = "OFF"

        self.gpib.write("SEL:CH%d: %s" % (channel, set))

    def set_impedance(self, channel, imped):
        """
        Sets the channel impedance
        :param channel: 1, 2, 3, or 4
        :param imped: 50 or 1e6
        :return:
        """

        if channel not in [1, 2, 3, 4]:
            print("Channel not correct. Doing nothing.")
            return

        if imped not in [50, 1e6]:
            print("Impedance option not correct. Doing nothing")
            return

        self.gpib.write("CH%d:TER %e" %(channel, imped))

    def set_vertical_range(self, channel, scale, offset=0):
        """
        Set the vertical range parameters for the specified channel. If they are
        None, nothing happens.
        :param channel: 1, 2, 3 or 4
        :param scale: for 50ohm impedance, 1x probe ratio: 1mV/div to 1 V/div; for 1Mohm impedance and 1x probe ratio: 1mV/div to 10V/div
        :param offset: see programming manual
        :return:
        """

        if channel not in [1, 2, 3 ,4]:
            print("Channel not correct. Doing nothing.")
            return

        if scale is not None:
            self.gpib.write("CH%d:SCA %.3f" % (channel, scale))

        if offset is not None:
            self.gpib.write("CH%d:OFFS %.3f" % (channel, offset))

    def set_horizontal_scale(self, scale):
        """
        Sets the horizontal scale in s/div
        :param scale: 200ps to 40s
        :return:
        """
        self.gpib.write("HOR:SCA %.7f" % scale)

    def read_horiz_rec_len(self):
        """
        Query horizontal record length
        :return: numpoints
        """
        numpoints = self.gpib.query("HOR:RECO?")
        # num_points = numpoints.strip() #for newer TDS7104
        num_points = str(numpoints[10:]).strip() #Works for TDS7104 with broken touchscreen, not newer one
        return float(num_points)

    def read_waveform(self, channels, file_name=None):
        """
        Reads the waveform on the screen in the specified channels
        :param channels: list with the channels whose waveform we w ant to obtain.
        :param file_name: if specified, it will save the data with the specified file name. Do not include the ".csv". One file for each channel
        :return: 2 lists, 2nd with n elemnts, where n is the number of channels.
                List 1: time
                List 2: [channel1_data, channel2_data, ...]
        """
        all_waveforms = []

        # Set to send ascii
        self.gpib.write("DATa:ENCdg ASCI")

        #Set number of points to be read to the horizontal record length (doesn't update automatically)
        horiz_rec_len = self.read_horiz_rec_len()
        self.gpib.write("DAT:STAR 1")
        self.gpib.write("DAT:STOP %d" %horiz_rec_len)

        for c in channels:
            if c not in [1, 2, 3, 4]:
                print("Specified channel not correct. Skipping it")
                continue

            # Choose source
            self.gpib.write("DATa:SOUrce CH%d" % c)

            self.gpib.write("CURV?")
            data = self.gpib.read_raw()

            wav_data_str = (str(data[6:-1].decode())).split(',')
            # print(wav_data_str)

            # wav_data = data
            wav_data = [float(i) for i in wav_data_str]
            read_points = len(wav_data)

            self.gpib.query("*ESR?")  # Clear standard event register
            err = self.gpib.query("EVMsg?")
            print(err)

            self.gpib.write("WFMOutpre?")
            preamble = self.gpib.read_raw()
            preamble_str = str(preamble[6:-1].decode()).split(';')
            # print(preamble_str)

            #Works for TDS7104 with broken touchscreen, not newer one
            num_points = float(preamble_str[6][5:])
            t_inc = float(preamble_str[9][4:])
            t_zero = float(preamble_str[10][4:])
            pt_offset = float(preamble_str[11][5:])
            v_mult = float(preamble_str[13][4:])
            v_offset = float(preamble_str[14][4:])
            v_zero = float(preamble_str[15][4:])

            #Create time array
            t = t_zero + (np.arange(num_points) - pt_offset)*t_inc
            print("Read %d points" % read_points)

            # Correct for offset and multiplying factor in data
            wav_data = list(v_zero + v_mult*(np.array(wav_data) - v_offset))

            # Save the data if necessary. Each channel will be stored in a different file
            if file_name is not None:

                # Create the csv file
                file_name_chan = file_name + "_channel_" + str(c) + ".csv"

                with open(file_name_chan, 'w+') as csvfile:
                    writer = csv.writer(csvfile)
                    writer.writerow(t)
                    writer.writerow(wav_data)

            all_waveforms.append(wav_data)

        return t, all_waveforms




if __name__ == '__main__':
    osc = TDS7104()
    osc.initialize()

    osc.read_waveform([1], 'unsync_OE1075.5_818-BB-35F_BOA-GSDev1a1253.3_Pulsed52ns9.5MHz')
    osc.close()

#%%
from photonmover.Interfaces.MSA import MSA
from photonmover.Interfaces.Instrument import Instrument

import pyvisa 
import time
import numpy as np
import matplotlib.pyplot as plt
import csv


GPIB_ADDRESS = "GPIB0::1::INSTR"


class YokogawaAQ6375(MSA, Instrument):
    def __init__(self):
        self.osa = None
    
    def initialize(self, GPIB_ADDRESS='GPIB0::1::INSTR', timeout=50000):
        """
        Connects to Yokogawa AQ6375 given resource_name (e.g., 'GPIB0::1::INSTR')
        """
        try:
            rm = pyvisa.ResourceManager()
            self.osa = rm.open_resource(GPIB_ADDRESS, timeout=timeout)
            print('Opening connnection to Yokogawa AQ6375 OSA')
        except BaseException:
            raise ValueError('Cannot connect to the Yokogawa AQ6375 OSA')
        # Verify connection
        print(f"Connected to: {self.osa.query('*IDN?').strip()}")
       
        
    def close(self):
        self.osa.close()
    
    def get_id(self):
        return ("YokogawaSpectrumAnalyzer")
            
    def read_data(self):
        """
        Performs a single sweep and returns wavelength and power data.
        tuple: (wavelengths, power_levels) as numpy arrays
        """
    
        # 1. Set to Single Sweep Mode and Trigger
        self.osa.write(":FORM:DATA ASCII")
        self.osa.write(":SENS:SWE:MODE SING") 
        self.osa.write("*CLS")                 # Clear status registers
        self.osa.write(":INIT:IMM")           # Trigger immediate sweep
        
        # 2. Wait for sweep completion
        # *OPC? returns '1' only after the sweep is finished
        print("Sweeping...")
        while True:
            try:
                if self.osa.query("*OPC?").strip() == "1":
                    break
            except pyvisa.errors.VisaIOError:
                # For very slow sweeps, the query might timeout; we keep waiting
                time.sleep(1)
        
        # 3. Retrieve Data
        # :TRACe:X? TRA returns wavelength data in meters
        # :TRACe:Y? TRA returns power data (usually dBm)
        # print("Transferring data...")
        raw_x = self.osa.query(":TRACe:X? TRA")
        raw_y = self.osa.query(":TRACe:Y? TRA")
        
        # print(f"raw_x: {raw_x}")
        raw_x = raw_x.replace("m","-") #Yokogawa randomly places m by some exponents
        # print(f"raw_y: {raw_y}")
        raw_y = raw_y.replace("/","-") # Y trace sometimes has '/' at begining
        
        
        # 4. Convert comma-separated strings to float arrays
        wavelengths = np.fromstring(raw_x, sep=',')
        powers = np.fromstring(raw_y, sep=',')
       
        # print(f"wavelengths: {wavelengths}")
        # print(f"powers: {powers}")
        
        return [wavelengths, powers]

    
    def set_rbw(self, rbw):
        self.osa.write(f":SENSe:BANDwidth:RESolution {rbw}NM")
        
    def set_center_wavelength(self, wavelength_nm):
        self.osa.write(f":DISP:TRACe:X:CENTer {wavelength_nm}NM")
    
    def set_span(self, span):
        self.osa.write(f":SENSe:WAV:SPAN {span}NM")
        
    def set_ref_level(self, ref_level):
        self.osa.write(f":DISP:TRACe:Y1:RLEV {ref_level}DBM")
        
    def set_num_points(self, num_points):
        self.osa.write(f"SENSe:SWEep:POINts {num_points}")
        
    def get_rbw(self):
        return self.osa.query_ascii_values(":SENS:BAND:RES?")[0]
    
    def get_center_wavelength(self):
        return self.osa.query_ascii_values(":DISP:TRACe:X:CENTer?")[0]
    
    def get_span(self):
        return self.osa.query_ascii_values(":SENSe:WAV:SPAN?")[0]
    
    def get_ref_level(self):
        return self.osa.query_ascii_values(":DISP:TRACe:Y1:RLEV?")[0]
    
    def get_num_points(self):
        return self.osa.query_ascii_values("SENSe:SWEep:POINts?")[0]
    
    def get_osa_parameters(self):
        """
        Read center wavelength, span, reference level, and RBW
        :return: [trace_len, osa_params] = dict(center_wavelength, span, reference_level, RBW, sensitivity)
        """
        center_wavelength = self.get_center_wavelength()
        span = self.get_span() #[m]
        reference_level = self.get_ref_level() #[dBm]
        RBW = self.get_rbw() #[m]
        # sensitivity = self.osa.query_ascii_values('SENS?') #[dBm]

        # osa_params = [center_wavelength, span, reference_level, RBW]
        osa_params = {
            'center wavelength':    center_wavelength,
            'span':                 span,
            'reference level':      reference_level,
            'rbw':                  RBW, 
            'sensitivity':          np.nan #For backward compatibility
        }

        trace_len = self.get_num_points()

        return [trace_len, osa_params]
        


#%%

if __name__ == '__main__':

    osa = YokogawaAQ6375()
    osa.initialize()
    [wavs, amps] = osa.read_data()
    print(wavs)
    print(amps)

    osa.close()


# %%

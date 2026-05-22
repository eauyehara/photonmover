import numpy as np
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from pint import UnitRegistry
from time import sleep
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, savgol_filter, hilbert, medfilt
from scipy.signal.windows import tukey
plt.rcParams['font.size'] = 14
u = UnitRegistry()
Q_ = u.Quantity

# Libraries for instrument control and saving data
from photonmover.instruments.Optical_spectrum_analyzers.HP70951B import HP70951B
from photonmover.instruments.Oscilloscopes.HP54750A import HP54750A
from photonmover.utils.hdf5_utils import *
from instrumental.drivers.daq.ni import Task # NIDAQ2,
from instrumental import instrument


data_dir = Path.home() / 'OneDrive\Documents\data\probe_station'
sample_dir = "BB_dispersion_setup"
logfile_name = "dispersion_log.csv"


#----------------------------------------------------------------------
# Placeholder for objective transmission - need to measure
T_obj = 1

# ds_noise = load_hdf5(data_dir / Path("supercont_osa_noiseFloor_rbw0.1nm_2026-04-23-01-20-51.h5"))

process_params = dict(
        cosφ_sgolay_window = 5,
        cosφ_polyorder = 3,
        num_periods = 3,
        medfilt_kernel = 5,
        good_fit = 0.8, 
        disp_polyorder = 3,            
    )

class dispersion_measurement:
    # code for controlling instruments for data acquisition
    
    def __init__(self):
        super().__init__()
        # Initialize class attributes here
        self.osa = None
        self.osc = None
        self.daq = None
        self.sample_shutter = None
        self.ref_shutter = None
    
    # -------------- Shutter control ------------------------------------------
    # Requires: daq
    # -------------------------------------------------------------------------
    def initialize_shutters(self):
        # Initialize daq
        self.daq = instrument("NIDAQ_USB-6259_2", reopen_policy='reuse')
        self.sample_shutter = self.daq.port1
        self.ref_shutter = self.daq.port0
        self.open_shutters()
        
    def open_shutters(self):
        self.ref_shutter.write(0x01)
        self.sample_shutter.write(0x01)

    def close_shutters(self):
        self.ref_shutter.write(0x00)
        self.sample_shutter.write(0x00)

    def sample_only(self):
        self.sample_shutter.write(0x01)
        self.ref_shutter.write(0x00)
        
    def ref_only(self):
        self.sample_shutter.write(0x00)
        self.ref_shutter.write(0x01)
    
 
     # -------------- MZI fringe measurement ------------------------------------------
     # Requires: daq (shutters), osa 
     # --------------------------------------------------------------------------------
    def initialize_osa(self):
        # Initialize osa
        self.osa = HP70951B()
        self.osa.initialize()
        
    def configure_osa(self, ref_level=Q_(-15, 'dBm'), rbw=0.1*u.nm, vbw=200*u.Hz, sensitivity=Q_(-50, 'dBm'), start_wl=750*u.nm, stop_wl=1450*u.nm, trace_length=2048, span=None):
        """
        Sets osa parameters. span not used, included to prevent errors in acquire spectra
        """
        self.osa.set_wl_axis(start_wl=f"{start_wl.to(u.nm).m} NM", end_wl=f"{stop_wl.to(u.nm).m} NM")
        self.osa.set_reference_level(ref_value=ref_level.m, ref_pos=None)
        self.osa.set_sensitivity(sens=sensitivity.m)
        self.osa.set_acq_bandwidth(res_bw = f"{rbw.to(u.nm).m} NM", video_bw = f"{vbw.to(u.Hz).m} HZ")
        self.osa.set_trace_length(trace_length=trace_length)
    
    
    def acquire_spectra(self, 
        osa_params = dict(
            start_wl = 750 * u.nm, 
            stop_wl = 1450 * u.nm, 
            span = 100 * u.nm, 
            rbw = 0.1 * u.nm, 
            vbw = 200 * u.Hz,
            sensitivity = Q_(-50, 'dBm'), 
            ref_level = Q_(-15, 'dBm'),
            trace_length = 2048,
        ), 
        fname=None, sample_dir=sample_dir, verbose=True, log = False):
        """
        Acquire a spectra by stitching together segments of length span.  
        Configures the osa with osa_params and records MZI spectra. Save to .h5 if fname is not None. If verbose=True, plots recorded spectra
        RETURN: ds of recorded MZI spectra
        """
        self.configure_osa(**osa_params)
        
        full_span = osa_params["stop_wl"] - osa_params["start_wl"]
        num_seg = int(full_span // osa_params["span"])
        
        # Create moving window of length span and store in array
        wavelength = np.zeros((num_seg * osa_params["trace_length"],))
        amplitude = np.zeros((num_seg * osa_params["trace_length"],))

        for i in range(num_seg):
            start_wl = osa_params["start_wl"] + i * osa_params["span"] 
            stop_wl = start_wl + osa_params["span"]
            self.osa.set_wl_axis(start_wl=f"{start_wl.to(u.nm).m} NM", end_wl=f"{stop_wl.to(u.nm).m} NM")

            print(f"start: {start_wl}")
            print(f"stop: {stop_wl} \n")
            [wavs, amps] = self.osa.read_data()
            
            wavelength[i*osa_params["trace_length"]: (i+1)*osa_params["trace_length"]] = wavs
            amplitude[i*osa_params["trace_length"]: (i+1)*osa_params["trace_length"]] = amps

        psd = amplitude - 10*np.log10(osa_params["rbw"].to(u.nm).m)

        if fname is not None:
            [_, meas_params] = self.osa.get_osa_parameters()

            time_stamp = new_timestamp()
            full_fname = fname + '_' + time_stamp + '.h5'
                
            fpath = data_dir / sample_dir / Path(full_fname)
            
            dump_hdf5(
                {
                "wavelength":   wavelength * u.m,
                "amplitude":    Q_(amplitude, 'dBm'),
                'psd':          psd, #dBm/nm
                "rbw":          meas_params["rbw"] * u.m,
                "sensitivity":   Q_(meas_params["sensitivity"], 'dBm'),
                "reference level":  Q_(meas_params['reference level'], 'dBm')
                }, 
                fpath=fpath, 
                open_mode='x'
            )
            
            if log:
                update_logfile(csv_logfile=logfile_name, fname=full_fname)
            
            ds = load_hdf5(fpath=fpath, verbose=False)
        else:
            ds = [wavelength, amplitude]
            
        if verbose:
            fig, ax = plt.subplots(1,1, figsize=(4,3), tight_layout=True)
            ax.plot(wavelength * 1e9, amplitude)
            ax.set_xlabel("Wavelength (nm)")
            ax.set_ylabel("Power (dBm)")
            
        return ds


    def acquire_normalized_spectra(self, osa_params, fname=None, sample_dir=sample_dir, verbose=True, log=True):
        """
        Sequentially acquires MZI spectra by stitching together segments of length span, sample arm spectra by closing reference shutter, and reference arm spectra by closing the sample shutter.
        If fname is not None, saves 3 files: MZI_fname..., SAMP_fname..., REF_fname...
        If verbose, plot:
            (1) Raw MZI, samp, and reference spectra vs wavelength 
            (2) Fringes normalized to supercontinuum spectrum
            (3) Sample Transmission
            (4) Extracted cosφ
        """
        # Create filenames
        mzi_fname = "MZI_" + fname
        samp_fname = "SAMP_" + fname
        ref_fname = "REF_" + fname
        
        # Acquire interference fringes - both shutters open
        self.open_shutters()
        ds_mzi = self.acquire_spectra(osa_params=osa_params, fname=mzi_fname, sample_dir=sample_dir, verbose=False, log=log)
        
        # Acquire sample arm only spectra - close reference shutter
        self.sample_only()
        ds_samp = self.acquire_spectra(osa_params=osa_params, fname=samp_fname, sample_dir=sample_dir, verbose=False, log=log)
        
        # Acquire ref arm only spectra - close sample shutter
        self.ref_only()
        ds_ref = self.acquire_spectra(osa_params=osa_params, fname=ref_fname, sample_dir=sample_dir, verbose=False, log=log)
        
        self.open_shutters()
        
        # Calculate normalized fringes
        amp_norm = dispersion_extractor.normalize_mzi(ds_mzi, ds_ref, ds_noise=None, verbose=False)
        
        # Calculate sample transmission
        T_samp = dispersion_extractor.calc_sample_T(ds_samp, ds_ref, ds_noise=None, verbose=False)
        
        # Extract cosφ
        ω_interp, amp_interp = dispersion_extractor.spec_to_freq(ds_mzi["wavelength"], amp_norm, nsamp=ds_mzi["wavelength"].shape[0], verbose=True)
        _, cosφ = dispersion_extractor.extract_cos(ω_interp, amp_interp, T_samp, verbose=False, smooth=False)
        
        λ_interp = (2*np.pi*u.c / ω_interp).to(u.nm)
        
        if verbose:
            fig, ax = plt.subplots(2,2, figsize=(8,6), tight_layout=True)
            ax[0,0].plot(ds_mzi["wavelength"].to(u.nm), ds_mzi["amplitude"], color='k')
            ax[0,0].plot(ds_samp["wavelength"].to(u.nm), ds_samp["amplitude"], color='b', label='S')
            ax[0,0].plot(ds_ref["wavelength"].to(u.nm), ds_ref["amplitude"], color='r', label='R')
            ax[0,0].set_xlabel("Wavelength (nm)")
            ax[0,0].set_ylabel("Power (dBm)")
            ax[0,0].legend()
            ax[0,0].set_title("Raw")
            ax[0,1].plot(λ_interp, amp_interp)
            ax[0,1].set_xlabel("Wavelength (nm)")
            ax[0,1].set_ylabel("Normalized (a.u.)")
            ax[0,1].set_title("Normalized Fringes")
            ax[1,0].plot(λ_interp, T_samp)
            ax[1,0].set_xlabel("Wavelength (nm)")
            ax[1,0].set_ylabel("Transmission")
            ax[1,0].set_title("Sample Transmission")
            ax[1,1].plot(λ_interp, cosφ)
            ax[1,1].set_xlabel("Wavelength (nm)")
            ax[1,1].set_ylabel("Normalized (a.u.)")
            ax[1,1].set_title("cosφ")
            
        return ds_mzi, ds_samp, ds_ref
    
    
    def measure_transmission(self, osa_params, fname=None, sample_dir=sample_dir, verbose=True):
        """
        Sequentially acquires device spectra by stitching together segments of length span, sample arm spectra by closing reference shutter, and reference arm spectra by closing the sample shutter.
        If fname is not None, saves 2 files: SAMP_fname..., REF_fname...
        If verbose, plot:
            (1) Raw samp and reference spectra vs wavelength 
            (2) Sample Transmission
        """
        # Create filenames
        samp_fname = "TSAMP_" + fname
        ref_fname = "TREF_" + fname
        
        # Acquire sample arm only spectra - close reference shutter
        self.sample_only()
        ds_samp = self.acquire_spectra(osa_params=osa_params, fname=samp_fname, sample_dir=sample_dir, verbose=False, log=False)
        
        # Acquire ref arm only spectra - close sample shutter
        self.ref_only()
        ds_ref = self.acquire_spectra(osa_params=osa_params, fname=ref_fname, sample_dir=sample_dir, verbose=False, log=False)
        
        # Calculate sample transmission (2x to account for second beamsplitter)
        T_samp = 2 * dispersion_extractor.calc_sample_T(ds_samp, ds_ref, ds_noise=None, verbose=False)
        
        
        if verbose:
            fig, ax = plt.subplots(2,1, figsize=(8,6), tight_layout=True)
            ax[0].plot(ds_samp["wavelength"].to(u.nm), ds_samp["amplitude"], color='b', label='S')
            ax[0].plot(ds_ref["wavelength"].to(u.nm), ds_ref["amplitude"], color='r', label='R')
            ax[0].set_xlabel("Wavelength (nm)")
            ax[0].set_ylabel("Power (dBm)")
            ax[0].legend()
            ax[0].set_title("Raw")
            ax[1].plot(ds_samp["wavelength"].to(u.nm), T_samp)
            ax[1].set_xlabel("Wavelength (nm)")
            ax[1].set_ylabel("Transmission")
            ax[1].set_title("Sample Transmission")
        
        self.open_shutters()

        return ds_samp, ds_ref
    
    
    @staticmethod
    def plot_spectra(fname_arr: list, sample_dir=sample_dir, label=None, ds_noise=None, xlim: tuple=None, ylim: tuple=None):
        fig, ax = plt.subplots(1, 1, figsize=(5,4), tight_layout=True)
        for ind, fname in enumerate(fname_arr): 
            fpath = data_dir / sample_dir / Path(fname)
            ds = load_hdf5(fpath=fpath)
            if label is not None:
                if len(label)==0:
                    ax.plot(ds['wavelength'].to(u.nm), ds['amplitude'].to(u.dBm), label=f'{ind}')
                else:
                    ax.plot(ds['wavelength'].to(u.nm), ds['amplitude'].to(u.dBm), label=label[ind])
                ax.legend(loc="upper right")
            else:
                ax.plot(ds['wavelength'].to(u.nm), ds['amplitude'].to(u.dBm))
            if ds_noise is not None:
                ax.plot(ds_noise['wavelength'].to(u.nm), ds_noise['amplitude'].to(u.dBm), color='silver')
            ax.set_ylabel('Power (dBm)')
            ax.set_xlabel("Wavelength (nm)")
            if xlim is not None:
                ax.set_xlim(xlim)
            if ylim is not None:
                ax.set_ylim(ylim)
        
        return fig
    
    # --------------- Sampling scope delay measurement ---------------------------------
    # Requires: daq (shutters), sampling scope
    # ----------------------------------------------------------------------------------
    def initialize_osc(self):
        self.osc = HP54750A()
        self.osc.initialize()
        

    def measure_delay(self, num_avg=16, wait=15*u.s, fname=None, sample_dir=sample_dir, sgolay_len=10, sgolay_poly=3, verbose=True, log=True):
        """
        Record sampling scope traces of the supercontinuum pulses in the reference arm and sample arms respectively for extraction of the sign of the group delay.
        If verbose, plot ref and sample traces on same time axis
        RETURN: delay between sample and reference peaks (φ1 = τ_samp - τ_ref)
        """
        # Configure osc num_avg
        self.osc.set_acq_averages(num_avg)
            
        # Sample delay
        self.sample_only()
        sleep(wait.m)
        t_samp, waveform_samp = self.osc.read_waveform([1]) 
        sleep(wait.m)
        t_samp = t_samp[0] * u.s
        waveform_samp = waveform_samp[0] * u.V
        
        # Reference delay
        self.ref_only()
        sleep(wait.m)
        t_ref, waveform_ref = self.osc.read_waveform([1]) 
        t_ref = t_ref[0] * u.s
        waveform_ref = waveform_ref[0] * u.V
        
        # if sgolay_len is not None:
        #     swaveform_samp = savgol_filter(waveform_samp, window_length=sgolay_len, polyorder=sgolay_poly)
        #     swaveform_ref = savgol_filter(waveform_ref, window_length=sgolay_len, polyorder=sgolay_poly)
        # else:
        #     swaveform_samp = waveform_samp
        #     swaveform_ref = waveform_ref
        
        # swaveform_ref *= u.V
        # swaveform_samp *= u.V
     
        # # calculate temporal delay between sample and reference pulses
        # ind_samp = np.argmax(swaveform_samp)
        # ind_ref = np.argmax(swaveform_ref)
        
        # τ_delay = t_samp[ind_samp] - t_ref[ind_ref]
        
        τ_delay = self.calc_delay(t_samp, t_ref, waveform_samp, waveform_ref, sgolay_len=sgolay_len, sgolay_poly=sgolay_poly, verbose=verbose)
        
        if fname is not None:
            time_stamp = new_timestamp()
            full_fname = 'DELAY_' + fname + '_' + time_stamp + '.h5'
            fpath = data_dir / sample_dir / Path(full_fname)
            
            dump_hdf5(
                {
                    "t_samp": t_samp,
                    "waveform_samp": waveform_samp,
                    "t_ref": t_ref,
                    "waveform_ref": waveform_ref,
                    "τ_delay": τ_delay
                }, 
                fpath=fpath, 
                open_mode='x'
            )
            if log:
                update_logfile(csv_logfile=logfile_name, fname=full_fname)
                
        #     # Shift to t_samp peak for plotting
        #     t0 = t_samp[ind_samp]
        #     t_samp_shift = (t_samp - t0).to(u.ps)
        #     t_ref_shift = (t_ref - t0).to(u.ps)
            
        #     # Normalize waveforms for plotting
        #     waveform_samp_norm = waveform_samp.to(u.mV) / np.max(waveform_samp.to(u.mV))
        #     waveform_ref_norm = waveform_ref.to(u.mV) / np.max(waveform_ref.to(u.mV))
        #     swaveform_samp_norm = swaveform_samp.to(u.mV) / np.max(swaveform_samp.to(u.mV))
        #     swaveform_ref_norm = swaveform_ref.to(u.mV) / np.max(swaveform_ref.to(u.mV))
      
        # if verbose:
        #     fig, ax = plt.subplots(1,1, figsize=(4,3), tight_layout=True)
        #     ax.scatter(t_samp_shift, waveform_samp_norm, color='b', label='S', marker='.')
        #     ax.plot(t_samp_shift, swaveform_samp_norm, color='b', alpha=0.8)
        #     ax.scatter(t_ref_shift, waveform_ref_norm, color='r', label='R', marker='.')
        #     ax.plot(t_ref_shift, swaveform_ref_norm, color='r', alpha=0.8)
        #     ax.vlines(x=[t_samp[ind_samp].to(u.ps).m, t_ref[ind_ref].to(u.ps).m] - t0.to(u.ps).m, ymin=np.min(waveform_ref_norm), ymax=1, linestyles=['--', '--'], colors=['k', 'k'])
        #     ax.set_xlabel("Time (ps)")
        #     ax.set_ylabel("Voltage (mV)")
        #     ax.set_title(f"τ_delay: {τ_delay.to(u.ps) :2.3f}")
            # ax.legend()
        
        self.open_shutters()
            
        return τ_delay
    
    
    @staticmethod
    def calc_delay(t_samp, t_ref, waveform_samp, waveform_ref, sgolay_len=10, sgolay_poly=3, verbose=True):
        """
        Plot sample and reference traces given ds_osc 
        """
        if sgolay_len is not None:
            swaveform_samp = savgol_filter(waveform_samp, window_length=sgolay_len, polyorder=sgolay_poly)
            swaveform_ref = savgol_filter(waveform_ref, window_length=sgolay_len, polyorder=sgolay_poly)
        else:
            swaveform_samp = waveform_samp
            swaveform_ref = waveform_ref
        
        swaveform_ref *= u.V
        swaveform_samp *= u.V
     
        # calculate temporal delay between sample and reference pulses
        ind_samp = np.argmax(swaveform_samp)
        ind_ref = np.argmax(swaveform_ref)
        
        τ_delay = t_samp[ind_samp] - t_ref[ind_ref]
        
        # Shift to t_samp peak for plotting
        t0 = t_samp[ind_samp]
        t_samp_shift = (t_samp - t0).to(u.ps)
        t_ref_shift = (t_ref - t0).to(u.ps)
        
        # Normalize waveforms for plotting
        waveform_samp_norm = waveform_samp.to(u.mV) / np.max(waveform_samp.to(u.mV))
        waveform_ref_norm = waveform_ref.to(u.mV) / np.max(waveform_ref.to(u.mV))
        swaveform_samp_norm = swaveform_samp.to(u.mV) / np.max(swaveform_samp.to(u.mV))
        swaveform_ref_norm = swaveform_ref.to(u.mV) / np.max(swaveform_ref.to(u.mV))
    
        if verbose:
            fig, ax = plt.subplots(1,1, figsize=(4,3), tight_layout=True)
            ax.scatter(t_samp_shift, waveform_samp_norm, color='b', label='S', marker='.')
            ax.plot(t_samp_shift, swaveform_samp_norm, color='b', alpha=0.8)
            ax.scatter(t_ref_shift, waveform_ref_norm, color='r', label='R', marker='.')
            ax.plot(t_ref_shift, swaveform_ref_norm, color='r', alpha=0.8)
            ax.vlines(x=[t_samp[ind_samp].to(u.ps).m, t_ref[ind_ref].to(u.ps).m] - t0.to(u.ps).m, ymin=np.min(waveform_ref_norm), ymax=1, linestyles=['--', '--'], colors=['k', 'k'])
            ax.set_xlabel("Time (ps)")
            ax.set_ylabel("Voltage (mV)")
            ax.set_title(f"τ_delay: {τ_delay.to(u.ps) :2.3f}")
            ax.legend()
            
        return τ_delay
        
    
# --------------- Utility Functions for creating/updating log files ---------------
def fname_to_log(fname):
    # Parse fname into dict for saving to logfile
    strs = fname.split('_')
    
    ftype = strs[0]
    print(ftype)
    if ftype == "DELAY":
        date = strs[4].replace('.h5','')
        params = ""
    else:
        stitch = strs[4]
        span = strs[5]
        rbw = strs[6]
        date = strs[7].replace('.h5','')
        params = stitch + '_' + span + '_' + rbw     
    chip = strs[1]
    w_dev = strs[2]
    L_dev = strs[3]

    log = dict(
        Plot = 0,
        ftype = ftype,
        Chip = chip,
        Length = L_dev,
        Width = w_dev,
        Params = params,
        Date = date,
        Notes = ""
    )
    return log


def update_logfile(fname, csv_logfile=logfile_name):
    # Create .csv logfile if doesn't exist.  Append info from fname to bottom
    try:
        df = pd.read_csv(csv_logfile)
        log = fname_to_log(fname)
        new_row = pd.DataFrame([log])
        df = pd.concat([df, new_row], ignore_index=True)
        df.to_csv(csv_logfile, index=False)
    except:
        log = fname_to_log(fname)
        df = pd.DataFrame([log])
        df.to_csv(csv_logfile, index=False)
      
        
#----------------------------------------------------------------------------------
class dispersion_extractor:
    # Code for extracting dispersion from dispersion measurement
    
    def __init__(self, ds_mzi, ds_samp, ds_ref, L_dev, devID, process_params: dict, ds_delay=None, ds_noise=None):
        super().__init__()
    
        # Normalize out supercontinuum intensity dependence from MZI output
        amp_norm = self.normalize_mzi(ds_mzi, ds_ref, ds_noise=ds_noise, verbose=False)
        N_interp = ds_mzi["wavelength"].shape[0]
        
        # Convert x axis to frequency and interpolate
        ω_interp, amp_interp = self.spec_to_freq(ds_mzi['wavelength'], amp_norm, nsamp=N_interp, verbose=True)
        λ_interp = (2*np.pi*u.c / ω_interp).to(u.nm)
        
        # Calculate sample transmission
        T_samp = self.calc_sample_T(ds_samp, ds_ref, ds_noise=ds_noise, verbose=False)
        
        # Convert sample transmission vs wavelength to transmission vs frequency
        _, T_samp_interp = self.spec_to_freq(ds_samp['wavelength'], T_samp, nsamp=N_interp, verbose=False)
        T = T_samp_interp * T_obj**2

        # Use calculated sample/objective transmission to extract cosφ term
        _, cosφ = self.extract_cos(ω_interp, amp_interp, T, verbose=True, smooth=True, window_length=5, polyorder=3)
        
        # If ds_delay is not None, extract the sign of the group delay difference
        if ds_delay is not None:
            φ1_sign = self.calc_φ1_sign(ds_delay=ds_delay)
        else:
            φ1_sign = None

        # Initialize class attributes here
        self.devID = devID              # string - description of dev
        self.L_dev = L_dev              # float - device length
        self.ω = ω_interp               # np.array - interpolated ω uniformly sampled in ω
        self.λ = λ_interp               # np.array - interpolated λ calculated from ω (non-uniform sampling)
        self.cosφ = cosφ                # np.array - extracted (smoothed) cosφ
        self.T_samp = T_samp_interp     # np.array - calculated sample transmission 
        self.norm_fringes = amp_interp  # np.array - calculated normalized fringes
        self.process_params = dict(                                     # nested dict - processing parameters 
            N_interp = N_interp,                                        # Number of points for interpolation of ω
            cosφ_sgolay_window = process_params["cosφ_sgolay_window"],  # Window size for sgolay filter of cosφ
            cosφ_sgolay_polyorder = process_params["cosφ_polyorder"],   # Polynomial order for sgolay filter of cosφ
            disp_polyorder = process_params["disp_polyorder"],          # Polynomial order for fitting phase before taking derivatives
            ω_lim = process_params["ω_lim"],                            # tuple - ω limits for phase fitting
            sliding_cos = dict(
                medfilt_kernel = process_params["medfilt_kernel"],      # for median filtering extracted phase before polynomial fitting
                num_periods = process_params["num_periods"],            # Number of periods to use in sliding cosine fit
                step = None,                                            # Step size for sliding cosine windows. If None, defaults to num_periods/4                              
                good_fit = process_params["good_fit"],                  # R^2 threshold of phase values to keep in sliding cosine fit
                ),
            hilbert = dict(
                tukey_alpha = process_params["tukey_alpha"],            # Apodization in hilbert transform
                trim_pts = process_params["trim_pts"]                           # Number of points to remove to remove hilbert transfrom edge artifacts
            )
        )  
        self.φ1_sign = φ1_sign          # string - '+' or '-' based on group delay difference measured with measure_delay()
        self.phase_results = dict(    # nested dict - phase extraction results for both ['sliding_cos'] and ['hilbert']
            sliding_cos = None,
            hilbert = None
            ) 
        self.GD = dict(                 # nested dict - calculated group delay
            sliding_cos = None, 
            hilbert = None
            )              
        self.GDD = dict(                # nested dict - calculated GDD 
            sliding_cos = None,
            hilbert = None
        )             
        self.GVD = dict(                # nested dict - calculated GVD
            sliding_cos = None,
            hilbert = None
        )             
        self.D = dict(                  # nested dict - calculated D 
            sliding_cos = None,
            hilbert = None
        )              
    
  # -------------- Methods for setting properties ------------------
  # Some of these properties are set upon creation of the dispersion_extractor object, but these functions can be used if a parameter needs to be modified
    def set_cosφ(self, smooth=True, window_length=10, polyorder=3, verbose=False, ω_lim=None): 
        self.ω, self.cosφ = self.extract_cos(ω=self.ω, norm_amp=self.norm_fringes, T=self.T_samp, smooth=smooth, window_length=window_length, polyorder=polyorder, verbose=verbose, ω_lim=ω_lim)
        self.λ = (2*np.pi*u.c / self.ω).to(u.nm)


    def set_process_params(self, process_params):
        for param, value in process_params.items():
            self.process_params[param] = value
          
          
    def set_φ1_sign(self, ds_delay=None, override=None):
        """
        Handle case where φ1_sign is zero?
        """
        if override is not None:
            if override == "-" or override == "+":
                self.φ1_sign = override
            else:
                raise ValueError("INVALID ENTRY! Specify '-' or '+'")
        elif ds_delay is not None:
            self.φ1_sign = self.calc_φ1_sign(ds_delay=ds_delay)
        else:
            raise ValueError("ERROR! Must provide either ds_delay or override value")    
        
            
    # -------------- Functions for correcting φ1 sign ------------------------------------
    @staticmethod
    def calc_φ1_sign(ds_delay):
        if ds_delay["τ_delay"] > 0:
            φ1_sign = "+"
        elif ds_delay["τ_delay"] < 0:
            φ1_sign = "-"
        else: # is 0
            φ1_sign = None
        return φ1_sign
    
    
    def fix_φ1_sign(self, φ1):
        """
        Cos symmetry creates ambiguity in sign of φ0 and φ1.  If self.φ1_sign is not None, enforce correct sign.
        RETURNS: -1 or 1
        """
    
        if self.φ1_sign is not None:
            if self.φ1_sign == "+" and  np.all(φ1[~np.isnan(φ1)] < 0):
                φ1_fix = -1  
            elif self.φ1_sign == "-" and np.all(φ1[~np.isnan(φ1)] > 0):
                φ1_fix = -1
            else:
                φ1_fix = 1
        else:
            φ1_fix = 1
        return φ1_fix
                
    # -------------- Functions for extracting cosφ ---------------------------------------
    # Modify to normalize out objective transmission xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    @staticmethod
    def calc_sample_T(ds_samp, ds_ref, ds_noise=None, verbose=False):
        """
        Calculate transmission through sample by dividing power in sample arm by power in ref arm 
        RETURN: T_lin (linear transmission), no units
        """
        # If wavelength axis is not same between ds_samp and ds_ref, raise error (Should be same from acquire_normalized_spectra(), so error in file loading if not)
        if np.all(ds_samp["wavelength"] != ds_ref["wavelength"]):
            raise ValueError("Wavelength axis of sample file doesn't match wavelength axis of reference file!")
        
        # Use psd to prevent errors from sample and ref acquisitions with different RBW (shouldn't happen with acquire_normalized_spectra())
        T_dB = ds_samp['psd'] - ds_ref['psd']
        T_lin = 10**(T_dB/10) #linear transmission
        
        if verbose:
            fig, ax = plt.subplots(2,1, figsize=(4,6), tight_layout=True)
            ax[0].plot(ds_ref['wavelength'].to(u.nm), ds_ref['psd'], color='k')
            ax[0].plot(ds_samp['wavelength'].to(u.nm), ds_samp['psd'], color='b')
            ax[1].plot(ds_samp['wavelength'].to(u.nm), T_dB, color='r')
        
            if ds_noise is not None:
                ax[0].plot(ds_noise['wavelength'].to(u.nm), ds_noise['psd'], color='silver')
                # ax[1].plot(ds_noise['wavelength'], ds_noise['psd']-ds_ref['psd'], color='silver')
            ax[0].set_ylabel("PSD (dBm/nm)")
            ax[1].set_xlabel("Wavelength (nm)")
            ax[1].set_ylabel("Transmission (dBm)")
            ax[0].set_xlim((np.min(ds_samp['wavelength'].to(u.nm)).m, np.max(ds_samp['wavelength'].to(u.nm)).m))
        return T_lin
    
    
    @staticmethod
    def normalize_mzi(ds_mzi, ds_ref, ds_noise=None, verbose=False):
        """
        Normalize MZI output intensity to intensity of supercontinuum (remove envelope on MZI fringes)
        RETURN: normalized LINEAR amplitude (no units)
        """
        # If wavelength axis is not same between ds_mzi and ds_ref, raise error (Should be same from acquire_normalized_spectra(), so error in file loading if not)
        if np.all(ds_mzi["wavelength"] != ds_ref["wavelength"]):
            raise ValueError("Wavelength axis of mzi file doesn't match wavelength axis of reference file!")
        
        amp_norm = ds_mzi['psd'] - ds_ref['psd']
        amp_norm = 10**(amp_norm/10)
        
        if verbose:
            fig, ax = plt.subplots(2,1, figsize=(4,6), tight_layout=True)
            ax[0].plot(ds_ref['wavelength'].to(u.nm), ds_ref['psd'], color='k')
            ax[0].plot(ds_mzi['wavelength'].to(u.nm), ds_mzi['psd'], color='b')
            ax[1].plot(ds_mzi['wavelength'].to(u.nm), amp_norm, color='r')
        
            if ds_noise is not None:
                ax[0].plot(ds_noise['wavelength'].to(u.nm), ds_noise['psd'], color='silver')
                # ax[1].plot(ds_noise['wavelength'], ds_noise['psd']-ds_ref['psd'], color='silver')
            ax[0].set_ylabel("PSD (dBm/nm)")
            ax[1].set_xlabel("Wavelength (nm)")
            ax[1].set_ylabel("Normalized")
            ax[0].set_xlim((np.min(ds_mzi['wavelength'].to(u.nm)).m, np.max(ds_mzi['wavelength'].to(u.nm)).m))
            ax[0].set_title("Normalized Fringes")
        return amp_norm
    
    
    @staticmethod
    def extract_cos(ω, norm_amp, T, verbose=False, smooth=True, window_length=10, polyorder=3, ω_lim: tuple =None, remove_baseline=True):
        """
        Extract cosine term of MZI output using analytic expression, supercontinuum-normalized MZI output, and total sample arm transmission. 
        Truncate range if ω_lim is not None.  Remove_baseline if remove_baseline=True.
        RETURN: (smoothed) cosφ
        """
        cosφ = 1 / (2*np.sqrt(T)) * (norm_amp - 1 - T)
        
        if smooth:
            cosφ_smooth = savgol_filter(cosφ, window_length=window_length, polyorder=polyorder)
            cos_out = cosφ_smooth
        else:
            cos_out = cosφ
            
        if ω_lim is not None:
            ω_trunc = ω[(ω > ω_lim[0]) & (ω < ω_lim[-1])]
            cosφ_trunc = cos_out[(ω > ω_lim[0]) & (ω < ω_lim[-1])]
        else:
            ω_trunc = ω
            cosφ_trunc = cos_out
            
        if remove_baseline:
            coeff = np.polyfit(ω_trunc.to(u.rad/u.fs).m, cosφ_trunc, deg=1)
            baseline = np.polyval(coeff, ω_trunc.to(u.rad/u.fs).m)
            cosφ_flat = cosφ_trunc - baseline
        else:
            cosφ_flat = cosφ_trunc
        
        if verbose:
            fig, ax = plt.subplots(1,1, figsize=(4,3), tight_layout=True)
            ax.plot(ω, cosφ, label='raw', color='r')
            if smooth:
                ax.plot(ω, cosφ_smooth, label='smooth', color='b')
            ax.plot(ω_trunc, cosφ_flat, color='m', alpha=0.7)
            ax.set_xlabel("ω (rad/s)")
            ax.set_ylabel("cosφ")
            ax.legend()
        return ω_trunc, cosφ_flat


    @staticmethod
    def spec_to_freq(wavelength, amplitude, nsamp=5000, verbose=True):
        """
        Convert wavelength axis to angular frequency and interpolate for uniform ω_samp
        Params: (Q_) wavelength, amplitude (Q_ or not)
        Return (Q_) ω_interp, amp_interp
        """
        
        # Handle amplitude units (if exists)
        try:
            amps = amplitude.m 
            un = amplitude.u
        except:
            amps = amplitude
            un = ''
        
        ω_raw = (2 * np.pi * u.c / (wavelength.to(u.nm).m * u.nm)).to(u.rad/ u.s)
        ω_interp = np.linspace(np.min(ω_raw), np.max(ω_raw), nsamp)
        amp_interp = np.interp(ω_interp.to(u.Hz).m, ω_raw[::-1].to(u.Hz).m, amps[::-1])
        
        if un != '': #add back units (if exists)
            amp_interp = Q_(amp_interp, amplitude.u)

        # Plot interpolated vs raw spectrum in freq
        if verbose:
            fig, ax = plt.subplots(1,1, figsize=(4,3), tight_layout=True)
            ax.plot(ω_raw, amps)
            ax.plot(ω_interp, amp_interp)
            ax.set_xlabel("ω (rad/s)")
            ax.set_ylabel(f"Amplitude ({un})")
            # ax.set_title("Freq Interpolated")
        return ω_interp, amp_interp


    # --------------- Sliding cosine functions ---------------------------------
    @staticmethod
    def mzi_func(ω, A, φ0, φ1): 
        φ = φ0 + φ1*ω 
        y = A * np.cos(φ)
        return y


    @staticmethod
    def get_Δω(ω, cosφ, verbose=False):
        """
        Find fringe spacing for initial fitting parameter guess.  Try max freq in fft first.  If fails
        Inputs: ω (Q_), amps (unitless)
        RETURN: Δω of fringe spacing, A amplitude of cos
        """
        try:
            if ω.shape[0] < 4:
                raise ValueError("Selected ω_segment has length < 4")

            dω = np.mean(np.diff(ω))

            fft_mag = np.abs(np.fft.rfft(cosφ)) / cosφ.shape[0] * 2
            freqs = np.fft.rfftfreq(n=ω.shape[0], d=dω)
            
            A = np.max(fft_mag) # Approximate magnitude of cos
            peak_ind = np.argmax(fft_mag[1:]) + 1 # skip DC
            Δω = 1/freqs[peak_ind]
            
            if verbose:
                fig, ax = plt.subplots(1,1, figsize=(4,3))
                ax.plot(1/freqs, fft_mag)
                ax.set_xlabel("ω (rad/s)")
                ax.set_ylabel("Power (W)")

        except:
            pk_ind, _ = find_peaks(cosφ, prominence=np.std(cosφ) * 0.3, distance=10, height=0.5*np.max(cosφ))
            if len(pk_ind) >= 2:
                Δω = np.median(np.diff(ω[pk_ind]))
                A = np.max(cosφ[pk_ind])
                print(f"Estimated Δω = {Δω:.3e}")
                
                if verbose:
                    fig, ax = plt.subplots(1,1, figsize=(4,3))
                    ax.plot(ω, cosφ)
                    ax.scatter(ω[pk_ind], cosφ[pk_ind], marker='x', color='k')
                    ax.set_xlabel("ω (rad/s)")
                    ax.set_ylabel("Power (W)")
            else:
                # Order of magnitude value for initial guess
                Δω = 1.00e-02 * u.rad/u.s
                print(f"find_peaks() can't extract fringe spacing")
        return Δω, A


    def extract_phase_sliding_cos(self, verbose=True):
        """
        Fit cosφ term, assuming polynomial phase.  
        PARAMS:
            - bounds: bounds for curve fit
            - window_len: sliding window (in rad/fs) used for fitting cosine
            - step: step size for sliding cosine window.  If none, defaults fo window_len / 4
            - good_fit: threshold of R² to include fit values of phase, else set to np.nan
        RETURN: dict of extracted phase coefficients.
        """
        
        if  self.φ1_sign is None:
            print("WARNING! The φ1 sign has not been set! The dispersion calculation may have a sign error!")
        
        # Calcualte the window length for the sliding cosine method given the specified number of periods
        pk_inds, _ = find_peaks(self.cosφ, prominence=np.std(self.cosφ) * 0.3, distance=10, height=0.5*np.max(self.cosφ))
        window_len = (self.ω[-1] - self.ω[0]).to(u.rad/u.fs).m / pk_inds.shape[0] * process_params["num_periods"]
        
        # Convert ω to rad/fs so large numbers don't break fitting
        ω_scale = (self.ω.to(u.rad/u.fs)).m #rad/fs
        good_fit = self.process_params['sliding_cos']['good_fit']
             
        if self.process_params["sliding_cos"]["step"] is None:
            step = window_len / 4
        
        ω_start = ω_scale[0]
        ω_end = ω_scale[-1]

        window_starts = np.arange(ω_start, ω_end - window_len + step * 0.1, step)
        n_windows = len(window_starts)

        if verbose:
            print(f"Number of windows: {n_windows}")
            print(f"Window length: {window_len:.3e} Hz")
            print(f"Step size: {step:.3e} Hz")   
            
        # Initialize arrays for fit variables
        ω_center = np.zeros(n_windows)
        amp = np.full(n_windows, np.nan)
        φ1 = np.full(n_windows, np.nan)
        φ0 = np.full(n_windows, np.nan)
        r_squared = np.zeros(n_windows)
        fit_success = np.zeros(n_windows, dtype=bool)

        # Full-length reconstructed fit (weighted average in overlaps)
        fit_full = np.zeros_like(self.cosφ)
        # fit_full = np.full(self.cosφ.shape, np.nan)
        weight_full = np.zeros_like(self.cosφ)
            
        for i, ws in enumerate(window_starts):
            we = ws + window_len

            # Extract window 
            mask = (ω_scale >= ws) & (ω_scale < we)
            ω_win = ω_scale[mask]
            cosφ_win = self.cosφ[mask]

            # Handle short segments at end of frequency vector (where we exceeds highest frequency)
            if len(ω_win) < 6:
                if verbose:
                    print(f"Window {i}: too few points ({len(ω_win)}), skipping")
                continue

            ω_center[i] = (ws + we) / 2.0
            
            # Initial phase needs to be close or fit will fail.
            Δω_win, A_win = self.get_Δω(ω=ω_win, cosφ=cosφ_win, verbose=False)
            φ1_est = (2*np.pi)/Δω_win     
            
            pk_inds, _ = find_peaks(cosφ_win, prominence=np.std(cosφ_win) * 0.3, distance=10, height=0.5*np.max(cosφ_win))
            if len(pk_inds) > 0:
                φ0_est = -φ1_est * ω_win[pk_inds[0]]
            else:
                φ0_est = 0.0
            p0 = [A_win, φ0_est, φ1_est] 

            bounds = (
                [0, -np.inf, 0.5*φ1_est],
                [np.inf, np.inf, 2*φ1_est]
            )   
            try:
                # Verify initial guess is finite
                y = self.mzi_func(ω_win, *p0)
                if not np.all(np.isfinite(y)):
                    raise ValueError("Non-finite initial residuals")
                
                popt, _ = curve_fit(self.mzi_func, ω_win, cosφ_win, p0=p0, method='trf', maxfev=100000, bounds=bounds)
            
                # R²
                y_fit = self.mzi_func(ω_win, *popt) 
                ss_res = np.sum((cosφ_win - y_fit)**2)
                ss_tot = np.sum((cosφ_win - np.mean(cosφ_win))**2)
                r_squared[i] = 1 - ss_res / ss_tot if ss_tot > 0 else 0
                fit_success[i] = True

                # Store results if r_squared > threshold
                if r_squared[i] > good_fit:
                    amp[i] = popt[0]
                    φ0[i] = popt[1]
                    φ1[i] = popt[2]
         
                
                # Accumulate into full fit (Hann window weighting for smooth stitching)
                hann = np.hanning(len(ω_win))
                fit_full[mask] += y_fit * hann
                weight_full[mask] += hann
                    
            except (RuntimeError, ValueError, TypeError) as e:
                if verbose:
                    print(f"  Window {i} at ω={ω_center[i]:.3e}: fit failed ({e})")
                fit_success[i] = False
            
        # Stitch full fit 
        valid = weight_full > 0
        fit_full[valid] /= weight_full[valid]
        fit_full[~valid] = np.nan
        
        # Enforce correct sign of φ1 if self.φ1_sign is not None
        if self.φ1_sign is not None:
            multiplier = self.fix_φ1_sign(φ1)
            φ1 *= multiplier
            φ0 *= multiplier    # If φ1 sign flipped, symmetric solution requires φ0 sign to be flipped
          
                
        # Store results in dict 
        results = {
            'method': 'sliding_cos',
            'ω_center': ω_center * u.rad/u.fs,
            'amp': amp,
            'φ0': φ0 * u.rad,
            'φ1': φ1 * u.fs,
            'r_squared': r_squared,
            'fit_success': fit_success,
            'cosφ_fit': fit_full,
        }
        
        self.phase_results["sliding_cos"] = results

        if verbose:
            n_ok = np.sum(fit_success)
            print(f"\nFit summary: {n_ok}/{n_windows} windows converged")
            print(f"Mean R²: {np.mean(r_squared[fit_success]):.4f}")
            ω_mean = np.mean(φ1[fit_success][~np.isnan(φ1)])
            print(f"Mean ω: {ω_mean:.6e} rad/s")
            
            fig = self.plot_sliding_cos()

        return results   


    def plot_sliding_cos(self):
        """
        Plot the sliding cosine fit results
        """
        results = self.phase_results['sliding_cos']
        good_fit_thresh = self.process_params['sliding_cos']['good_fit']
        mask = results['fit_success']

        fig, ax = plt.subplots(3, 1, figsize=(6, 10), sharex=False, tight_layout=True)
        # Raw data vs reconstructed fit
        ax[0].plot(self.ω.to(u.rad/u.fs), self.cosφ, 'b-', alpha=0.7, label='raw')
        ax[0].plot(self.ω.to(u.rad/u.fs), results['cosφ_fit'], 'r-', alpha=0.8, label='fit')
        ax[0].set_ylabel('Normalized (a.u.)')
        ax[0].set_xlabel('ω (rad/fs)')
        ax[0].legend()
        ax[0].set_title('Sliding cosφ fit')
        # Extracted Instantaneous frequency φ1
        ax[1].plot(results['ω_center'][mask].to(u.rad/u.fs), results['φ1'][mask],'go-', ms=3)
        ax[1].set_ylabel('φ1 (fs)')
        ax[1].set_xlabel('ω (rad/fs)')
        ax[1].set_title('Instantaneous Frequency')
        # R² Fit quality
        ax[2].plot(results['ω_center'][mask].to(u.rad/u.fs), results['r_squared'][mask],'ms-', ms=2)
        ax[2].plot(results['ω_center'][mask].to(u.rad/u.fs), good_fit_thresh*np.ones(results['ω_center'][mask].shape[0]), color='k', linestyle='--')
        ax[2].set_ylabel('R²')
        ax[2].set_xlabel('ω (rad/fs)')
        ax[2].set_title('Fit Quality per Window')
        ax[2].set_ylim([0, 1.05])
        # plt.show()
        return fig


    # ------------- Hilbert Transform functions -------------------------------
    def extract_phase_hilbert(self, verbose=True):
        """
        Extract phase φ(ω) from MZI fringes using the Hilbert transform.
        RETURN: result (dict)
                'φ'              : (Q_) unwrapped spectral phase φ(ω)
                'amp'            : instantaneous amplitude
                'analytic_signal': (Q_) complex analytic signal 
                'trim_slice'     : slice for removing Hilbret artifacts       
        """
        if  self.φ1_sign is None:
            print("WARNING! The φ1 sign has not been set! The dispersion calculation may have a sign error!")
            
        ω_scale = self.ω.to(u.rad/u.fs).m
        params = self.process_params['hilbert']

        # Hilbert transform to create complex signal 
        w = tukey(len(self.cosφ), alpha=params['tukey_alpha'])
        analytic = hilbert(self.cosφ * w)
        amp_inst = np.abs(analytic)
        φ_inst = np.unwrap(np.angle(analytic))

        if verbose:
            fig, ax = plt.subplots(1,1, figsize=(4,3))
            ax.plot(ω_scale, np.real(analytic), label='R', color='blue')
            ax.plot(ω_scale, np.imag(analytic), label='I', color='red')
            plt.legend()

        # Trim edges (Hilbert artifacts)
        # trim = max(len(ω_scale) // 50, params["trim_pts"])
        trim = params["trim_pts"]
        sl = slice(trim, -trim)

        results = {
            'φ':                φ_inst * u.rad,
            'amp':              amp_inst,
            'analytic_signal':  analytic,
            'trim_slice':       sl,
        }
        
        self.phase_results["hilbert"] = results

        return results


    # ------------ Dispersion Calculation Functions ----------------
    def calc_dispersion(self, verbose=True):
        """
        Calculate GDD (or GVD if L_dev is not None) and D.
        If self.phase_results is None, calls extract_phase_sliding_cos() and extract_phase_hilbert().
        Calculates dispersion for both sliding cos and hilbert transform methods.
        Sets self.GVD, self.GDD, and self.D for (each) method
        """
        sc_params = self.process_params["sliding_cos"]
        disp_polyorder = self.process_params["disp_polyorder"]
        trim = self.process_params["hilbert"]["trim_pts"]
        
        # ---------------- Check if self.phase_results has data --------------------
        if self.phase_results["sliding_cos"] is None:
            sc_results = self.extract_phase_sliding_cos()
        elif self.phase_results["hilbert"] is None:
            h_results = self.extract_phase_hilbert()
        else: # Already solved for phase
            sc_results = self.phase_results["sliding_cos"]
            h_results = self.phase_results["hilbert"]
        
        # ---------------- Sliding cosine dispersion extraction ----------------------
        # Construct φ(ω) and its derivatives
        φ1_sc = sc_results["φ1"].to(u.fs).m
        mask = np.isfinite(φ1_sc) #remove nan from fit

        # Interpolate back to ω_center after removing nans
        dφdω_sc = np.interp(sc_results["ω_center"].to(u.rad/u.fs).m, sc_results["ω_center"].to(u.rad/u.fs).m[mask], φ1_sc[mask])
        
        # Median filter 
        if sc_params["medfilt_kernel"] is not None:
            dφdω_filt = medfilt(dφdω_sc, kernel_size=sc_params["medfilt_kernel"])
        else:
            dφdω_filt = dφdω_sc
        
        # Fit to polynomial
        φ1_coeffs_sc = np.polyfit(sc_results["ω_center"].to(u.rad/u.fs).m, dφdω_filt, deg=disp_polyorder-1) # 1 order polynomial order than Hilbert transform since fitting φ1, not φ
        self.GD["sliding_cos"] = np.polyval(φ1_coeffs_sc, self.ω.to(u.rad/u.fs).m) * u.rad/u.fs
        
        # Construct derivatives
        φ2_coeffs_sc = np.polyder(φ1_coeffs_sc, 1) # φ2
        self.GDD["sliding_cos"] = np.polyval(φ2_coeffs_sc, self.ω.to(u.rad/u.fs).m) * u.fs**2
        
        # ----------------- Hilbert Transform dispersion extraction --------------------------
        sl = h_results['trim_slice']
        
        dω = np.mean(np.diff(self.ω.to(u.rad/u.fs).m))
        
        # Remove hilbert transform edge artifacts 
        coeffs_h = np.polyfit(self.ω[sl].to(u.rad/u.fs).m, h_results["φ"][sl].to(u.rad).m, deg=disp_polyorder) 
        φ_poly_h = np.polyval(coeffs_h, self.ω[sl].to(u.rad/u.fs).m)
    
        # Calculate derivatives
        gd_h = np.gradient(φ_poly_h, dω)                 # dφ/dω
        gdd_h = np.gradient(gd_h, dω)                    # d²φ/dω²
        
        # Enforce correct sign of group delay
        multiplier = self.fix_φ1_sign(gd_h)
        self.phase_results['hilbert']['φ'] *= multiplier
        
        # Pad with Nans to recover full size of ω
        self.GD["hilbert"] = np.concatenate((np.full(φ_poly_h[0:trim].shape, np.nan), multiplier*gd_h, np.full(φ_poly_h[-trim:].shape, np.nan))) * u.rad/u.fs    
        self.GDD["hilbert"] = np.concatenate((np.full(φ_poly_h[0:trim].shape, np.nan), multiplier*gdd_h, np.full(φ_poly_h[-trim:].shape, np.nan))) * u.fs**2  
        
        
        if self.L_dev is not None:
            # GVD
            self.GVD["sliding_cos"] = (1/self.L_dev * self.GDD["sliding_cos"]).to(u.fs**2 / u.mm)
            self.GVD["hilbert"] = (1/self.L_dev * self.GDD["hilbert"]).to(u.fs**2 / u.mm)
            y1_sc = self.GVD["sliding_cos"]
            y1_h = self.GVD["hilbert"]
            y1_label = 'GVD (fs²/mm)'
            
            # D parameter
            self.D["sliding_cos"] = (-self.GVD["sliding_cos"]* 2*np.pi*u.c / self.λ**2).to(u.ps / (u.nm * u.km))
            self.D["hilbert"] = (-self.GVD["hilbert"]* 2*np.pi*u.c / self.λ**2).to(u.ps / (u.nm * u.km))
            y2_sc = self.D["sliding_cos"]
            y2_h = self.D["hilbert"]
            y2_label = 'D (ps/(nm km))'
        else:
            # GDD
            y1_sc = self.GDD["sliding_cos"].to(u.fs**2)
            y1_h = self.GDD["hilbert"].to(u.fs**2)
            y1_label = 'GDD (fs²)'
            
            # D*L_dev
            y2_sc = (-y1_sc * 2*np.pi*u.c / self.λ**2).to(u.ps / u.nm)
            y2_h = (-y1_h * 2*np.pi*u.c / self.λ**2).to(u.ps / u.nm)
            y2_label = 'D (ps/nm)'

        # Print zdw
        zdw_ind_sc = np.nanargmin(np.abs(y1_sc))
        zdw_ind_h = np.nanargmin(np.abs(y1_h))
        zdw_sc = self.λ[zdw_ind_sc].to(u.nm)
        zdw_h = self.λ[zdw_ind_h].to(u.nm)
        print(f"ZDW (sliding cos) = {zdw_sc :2.3f}")
        print(f"ZDW (hilbert) = {zdw_h :2.3f}")
        
        if verbose:
            fig, ax = plt.subplots(3,1, figsize=(4,6), tight_layout=True)
            ax[0].scatter(sc_results["ω_center"].to(u.rad/u.fs), dφdω_sc, marker='.', color='k')
            ax[0].scatter(sc_results["ω_center"].to(u.rad/u.fs), dφdω_filt, marker='.', color='r')
            ax[0].plot(self.ω.to(u.rad/u.fs), self.GD["sliding_cos"], color='b', label='sc', alpha=0.8)
            ax[0].plot(self.ω.to(u.rad/u.fs), self.GD["hilbert"], color='m', label='h', alpha=0.7)
            ax[0].legend()
            ax[0].set_xlabel("ω (rad/fs)")
            ax[0].set_ylabel("dφ/dω (fs)")
            
            ax[1].plot(self.ω.to(u.rad/u.fs), y1_sc, color='b', label='sc', alpha=0.8)
            ax[1].plot(self.ω.to(u.rad/u.fs), y1_h, color='m', label='h', alpha=0.7)
            ax[1].plot(self.ω.to(u.rad/u.fs), np.zeros(self.ω.shape), linestyle='--', color='k')
            ax[1].legend()
            ax[1].set_ylabel(y1_label)
            ax[1].set_xlabel("ω (rad/fs)")
            
            ax[2].plot(self.λ.to(u.nm), y2_sc, color='b', label='sc', alpha=0.8)
            ax[2].plot(self.λ.to(u.nm), y2_h, color='m', label='h', alpha=0.7)
            ax[2].plot(self.λ.to(u.nm), np.zeros(self.λ.shape), linestyle='--', color='k')
            ax[2].legend()
            ax[2].set_xlabel("λ (nm)")
            ax[2].set_ylabel(y2_label)
            
        return fig, ax


    def plot_dispersion(self, fig=None, ax=None, wavelength_axis=True):
        """
        Plot the extracted dispersion for both methods
        """
        if wavelength_axis:
            x_axis = self.λ.to(u.nm)
            xlabel = "Wavelength (nm)"
        else:
            x_axis = self.ω.to(u.rad/u.fs)
            xlabel = "ω (rad/fs)"
            
        if self.GVD is None:
            if fig is None and ax is None:
                fig, ax = plt.subplots(1,1, figsize=(4,3), tight_layout=True)
            
            ax.plot(x_axis, self.GDD["sliding_cos"], color='b', label='sc')  
            ax.plot(x_axis, self.GDD["hilbert"], color='r', label='h') 
            ax.plot(self.ω.to(u.rad/u.fs), np.zeros(self.ω.shape), linestyle='--', color='k')
            ax.legend()
            ax.set_ylabel("GDD (fs²)")
            ax.set_xlabel(xlabel)
        else:
            if fig is None and ax is None:
                fig, ax = plt.subplots(2,1, figsize=(4,3), tight_layout=True)
                
            ax[0].plot(x_axis, self.GVD["sliding_cos"], color='b', label='sc')  
            ax[0].plot(x_axis, self.GVD["hilbert"], color='r', label='h') 
            ax[0].plot(x_axis, np.zeros(self.ω.shape), linestyle='--', color='k')
            ax[0].legend()
            ax[0].set_ylabel("GVD (fs²/mm)")
            ax[0].set_xlabel(xlabel)
            
            ax[1].plot(x_axis, self.D["sliding_cos"], color='b', label='sc')  
            ax[1].plot(x_axis, self.D["hilbert"], color='r', label='h') 
            ax[1].plot(x_axis, np.zeros(self.ω.shape), linestyle='--', color='k')
            ax[1].legend()
            ax[1].set_ylabel("D (ps / (nm km))")
            ax[1].set_xlabel(xlabel)

        return fig, ax
        
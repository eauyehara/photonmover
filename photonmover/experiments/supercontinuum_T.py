#%%
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from pint import UnitRegistry
from photonmover.utils.hdf5_utils import *
from photonmover.instruments.Optical_spectrum_analyzers.HP70951B import HP70951B

u = UnitRegistry()
Q_ = u.Quantity

osa = HP70951B()
osa.initialize()

# %% Configure osa
ref_level = Q_(-35, 'dBm') #[dBm]
rbw = 1 * u.nm #[nm]
sensitivity = Q_(-50, 'dBm') #[dBm]
span = None #800 * u.nm 
center = None #1300 * u.nm 
start_wl = 600 * u.nm 
end_wl = 1700 * u.nm 

osa.set_wl_axis(start_wl=f"{start_wl.to(u.nm).m} NM", end_wl=f"{end_wl.to(u.nm).m} NM")
osa.set_reference_level(ref_value=ref_level.m, ref_pos=None)
osa.set_sensitivity(sens=sensitivity.m)
osa.set_acq_bandwidth(res_bw = f"{rbw.to(u.nm).m} NM")

#%% Record OSA NOISE FLOOR trace and osa parameters.  Save to hdf5 >>>>>>>>>>>>>>>>>>>

fname = "supercont_osa_noiseFloor"

data_dir = Path.home() / 'OneDrive\Documents\data\probe_station'
time_stamp = new_timestamp()
fpath = data_dir / Path(fname + '_' + time_stamp + '.h5')

[wavs1, amps1] = osa.read_data()
[_, osa_params] = osa.get_osa_parameters()

psd1 = amps1 - 10*np.log10(rbw.to(u.nm).m)

dump_hdf5({
    "wavelength":   wavs1 * u.m,
    "amplitude":    Q_(amps1, 'dBm'),
    'psd':          psd1,
    "rbw":          osa_params["rbw"] * u.m,
    "sensitivity":   Q_(osa_params["sensitivity"], 'dBm'),
    "reference level":  Q_(osa_params['reference level'], 'dBm')
    }, 
    fpath=fpath, 
    open_mode='x'
)

#%% Record PROBE-to-PROBE trace and osa parameters.  Save to hdf5 >>>>>>>>>>>>>>>>>>>

fname = "supercont_probe-to-probe_3"
# fname = "supercont_atten"

data_dir = Path.home() / 'OneDrive\Documents\data\probe_station'
time_stamp = new_timestamp()
fpath = data_dir / Path(fname + '_' + time_stamp + '.h5')

[wavs0, amps0] = osa.read_data()
[_, osa_params] = osa.get_osa_parameters()

psd0 = amps0 - 10*np.log10(rbw.to(u.nm).m)

dump_hdf5({
    "wavelength":   wavs0 * u.m,
    "amplitude":    Q_(amps0, 'dBm'),
    'psd':          psd0,
    "rbw":          osa_params["rbw"] * u.m,
    "sensitivity":   Q_(osa_params["sensitivity"], 'dBm'),
    "reference level":  Q_(osa_params['reference level'], 'dBm')
    }, 
    fpath=fpath, 
    open_mode='x'
)

#%% Record DEVICE trace and osa parameters. Save to hdf5 >>>>>>>>>>>>>>>>>>>>>
fname = "supercont_Wafer8_R2C3_W2.0um_ref"
# fname = "test"

data_dir = Path.home() / 'OneDrive\Documents\data\probe_station'
time_stamp = new_timestamp()
fpath = data_dir / Path(fname + '_' + time_stamp + '.h5')

[wavs, amps] = osa.read_data()
[_, osa_params] = osa.get_osa_parameters()

psd = amps - 10*np.log10(rbw.to(u.nm).m)

dump_hdf5({
    "wavelength":   wavs * u.m,
    "amplitude":    Q_(amps, 'dBm'),
    'psd':          psd,
    "rbw":          osa_params["rbw"] * u.m,
    "sensitivity":   Q_(osa_params["sensitivity"], 'dBm'),
    "reference level":  Q_(osa_params['reference level'], 'dBm')
    }, 
    fpath=fpath, 
    open_mode='x'
)

#%% Load probe-to-probe from file
fname = "supercont_probe-to-probe_2025-10-21-14-26-53.h5"
fpath = data_dir / Path(fname)
ds = load_hdf5(fpath=fpath)

# fname2 = "supercont_atten_2025-10-17-18-49-10.h5"
# fpath2 = data_dir / Path(fname2)
# ds2 = load_hdf5(fpath=fpath2)

# amps=ds2["amplitude"].m
# wavs=ds2["wavelength"].m
#%% Plot device Transmission

#Assuming wavelength axis same for probe-to-probe and device measurments
T = np.array(amps) - ds['amplitude'].m

fig, ax = plt.subplots(3, 1, figsize=(4,8), tight_layout=True)
ax[0].plot(wavs*1e9, amps)
ax[0].plot(wavs1*1e9, amps1, color='silver')
ax[0].set_ylabel('Device PSD (dBm/nm)')
ax[1].plot(ds['wavelength'].to(u.nm), ds['amplitude'])
ax[1].set_ylabel('Probe PSD (dBm/nm)')
ax[2].plot(wavs*1e9, T)
ax[2].set_ylabel('Transmission (dB)')
# ax[2].set_ylim((-25,0.5))
ax[2].set_xlabel("Wavelength (nm)")

#%% Save Transmission plot >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

fig_name = data_dir / Path(fname + '_' + time_stamp + '.png')

fig, ax = plt.subplots(1, 1, figsize=(4,3), tight_layout=True)
ax.plot(wav, T)
ax.set_ylabel('Transmission (dB)')
ax.set_xlabel("Wavelength (nm)")

#%% Post-processing plotting >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# from scipy.signal import savgol_filter

fname = "supercont_probe-to-probe_2025-10-21-14-26-53.h5"
fname0 = "supercont_probe-to-probe_2_2025-10-21-14-34-20.h5"
fname1 = "supercont_probe-to-probe_3_2025-10-21-14-36-49.h5"
fpath = data_dir / Path(fname)
fpath0 = data_dir / Path(fname0)
fpath1 = data_dir / Path(fname1)
ds0 = load_hdf5(fpath=fpath)
ds1 = load_hdf5(fpath=fpath0)
ds2 = load_hdf5(fpath=fpath1)

plt.figure(figsize = (6,4), tight_layout=True)
plt.plot(ds0["wavelength"].to(u.nm).m, ds0["amplitude"])
plt.plot(ds1["wavelength"].to(u.nm).m, ds1["amplitude"])
plt.plot(ds2["wavelength"].to(u.nm).m, ds2["amplitude"])
plt.xlabel("Wavelength (nm)")
plt.ylabel("Power (dBm)")

plt.figure(figsize = (6,4), tight_layout=True)
window_length = 20
polyorder = 5
plt.plot(ds0["wavelength"].to(u.nm).m, savgol_filter(ds0["amplitude"], window_length=window_length, polyorder=polyorder))
plt.plot(ds1["wavelength"].to(u.nm).m, savgol_filter(ds1["amplitude"], window_length=window_length, polyorder=polyorder))
plt.plot(ds2["wavelength"].to(u.nm).m, savgol_filter(ds2["amplitude"], window_length=window_length, polyorder=polyorder))
plt.xlabel("Wavelength (nm)")
plt.ylabel("Power (dBm)")

#%%
p2p_smooth = savgol_filter(ds2["amplitude"], window_length=window_length, polyorder=polyorder)
amps_smooth = savgol_filter(amps, window_length=window_length, polyorder=polyorder)
T = savgol_filter(np.array(amps), window_length=window_length, polyorder=polyorder) - p2p_smooth

fig, ax = plt.subplots(3, 1, figsize=(4,8), tight_layout=True)
ax[0].plot(wavs*1e9, amps_smooth)
ax[0].plot(wavs1*1e9, amps1, color='silver')
ax[0].set_ylabel('Device PSD (dBm/nm)')
ax[1].plot(ds['wavelength'].to(u.nm), p2p_smooth)
ax[1].set_ylabel('Probe PSD (dBm/nm)')
ax[2].plot(wavs*1e9, T)
ax[2].set_ylabel('Transmission (dB)')
# ax[2].set_ylim((-25,0.5))
ax[2].set_xlabel("Wavelength (nm)")
# %%

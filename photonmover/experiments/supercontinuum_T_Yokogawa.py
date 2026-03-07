#%%
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
# from pint import UnitRegistry
from photonmover.utils.hdf5_utils import *
from photonmover.instruments.Optical_spectrum_analyzers.YokogawaAQ6375 import YokogawaAQ6375

u = UnitRegistry()
Q_ = u.Quantity

osa = YokogawaAQ6375()
osa.initialize()
# %%
ref_level = Q_(-35, 'dBm') #[dBm]
rbw = 0.1 * u.nm #[nm]
vbw = 100 * u.Hz
sensitivity = Q_(-50, 'dBm') #[dBm]
span = None #800 * u.nm 
center = None #1300 * u.nm 
start_wl = 1550 * u.nm 
end_wl = 1800 * u.nm 


osa.set_rbw(rbw.to(u.nm).m)
osa.set_ref_level(ref_level.m)
span = end_wl - start_wl
center_wl = np.round(span/2) + start_wl
osa.set_span(span.to(u.nm).m)
osa.set_center_wavelength(center_wl.to(u.nm).m)
#%% Record OSA NOISE FLOOR trace and osa parameters.  Save to hdf5 >>>>>>>>>>>>>>>>>>>

fname = "supercont_Yokogawa_noiseFloor"

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


#%% Record DEVICE trace and osa parameters. Save to hdf5 >>>>>>>>>>>>>>>>>>>>>
fname = "Yokogawa_supercont_Wafer8_R2C3_Ring8_W0.8um_gap0.55um_rbw0.1nm_1550-1800nm_1mW"
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
fig, ax = plt.subplots(1, 1, figsize=(6,4), tight_layout=True)
ax.plot(wavs*1e9, amps)
ax.plot(wavs1*1e9, amps1, color='silver')
ax.set_ylabel('Power(dBm)')
ax.set_xlim((min(wavs)*1e9, max(wavs)*1e9))
ax.set_ylim((-100, -60))
# %%

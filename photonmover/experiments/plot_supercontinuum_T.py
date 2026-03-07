"""
Load and plot supercontinuum Transmission files
"""
#%%
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from pint import UnitRegistry
from photonmover.utils.hdf5_utils import *
u = UnitRegistry()
Q_ = u.Quantity

#%% Load probe-to-probe from file >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
data_dir = Path.home() / 'OneDrive\Documents\data\probe_station' 


# 3/3/26 Tarik Alluxa - no objectives
sample_dir = "Tarik_Alluxa"
fname_noise = ["supercont_osa_noiseFloor_2026-03-03-14-37-22.h5"]
fname_device = ["Alluxa895.8nm_BPF_600-1700nm_2026-03-03-14-42-18.h5",
                "Alluxa830nm_LPF_600-1700nm_2026-03-03-14-58-44.h5"]
fname_probe = ["supercont_probe-to-probe_noOBJ_2026-03-03-14-39-19.h5"]
# dev_label = ["SBS15A_1550"]
dev_label = ["BPF", "LPF"]

# ----------------------------------------------------------------------
# # 2/23/26 Dilute SiN SBS15A - 0.95 NA objectives
# sample_dir = "LL_dilute_SiN"
# fname_noise = ["supercont_osa_noiseFloor_2026-02-22-23-32-03.h5"]
# fname_device = ["SBS15A_1550_rbw1nm_600-1700nm_vbw400Hz_2026-02-22-23-30-10.h5",
#                 "SBS15A_1550_AC050-008-C_rbw1nm_600-1700nm_vbw400Hz_2026-02-23-20-59-47.h5",
#                 "SBS15A_1550_AC050-008-C_FESH0950_rbw1nm_600-1700nm_vbw400Hz_2026-02-23-21-07-43.h5",
#                 "SBS15A_1550_AC050-008-C_SiPM_rbw1nm_600-1700nm_vbw400Hz_2026-02-23-21-16-51.h5",
#                 "SBS15A_1550_SiPM_rbw1nm_600-1700nm_vbw400Hz_2026-02-22-23-35-52.h5"]
# fname_probe = ["supercont_probe-to-probe_0_2026-02-22-22-10-08.h5",
#                "supercont_probe-to-probe_AC050-008-C_2026-02-23-20-26-57.h5"]
# # dev_label = ["SBS15A_1550"]
# dev_label = ["OBJ0.95NA", "ACH", "ACH_SP950", "ACH_SiPM", "OBJ0.95NA_SiPM"]

# -----------------------------------------------------------------------
# 2/15/26 Wafer8 R2C3 Rings - 0.95 NA objectives
# sample_dir = 'LL_SiN_Rings'
# fname_noise = ["supercont_osa_noiseFloor_2026-02-15-17-36-43.h5", # 600-1700 nm
#                "supercont_osa_noiseFloor_2026-02-15-22-40-02.h5"] # 1550-1700 nm
# fname_superc = ["supercont_rbw1nm_600-1700nm_vbw100Hz_2026-02-16-02-02-22.h5"]

# fname_device_hires_HP = [
#                 "supercont_Wafer8_R2C3_Ring2_W0.8um_gap0.3um_rbw0.1nm_1550-1700nm_vbw100Hz_1mW_2026-02-15-23-43-21.h5",
#                 "supercont_Wafer8_R2C3_Ring3_W0.8um_gap0.35um_rbw0.1nm_1550-1700nm_vbw100Hz_1mW_2026-02-16-00-04-20.h5",
#                 "supercont_Wafer8_R2C3_Ring4_W0.8um_gap0.4um_rbw0.1nm_1550-1700nm_vbw100Hz_1mW_2026-02-16-00-36-24.h5",
#                 "supercont_Wafer8_R2C3_Ring5_W0.8um_gap0.45um_rbw0.1nm_1550-1700nm_vbw100Hz_1mW_2026-02-16-00-56-04.h5",
#                 "supercont_Wafer8_R2C3_Ring7_W0.8um_gap0.55um_rbw0.1nm_1550-1700nm_vbw100Hz_1mW_2026-02-16-01-30-36.h5",
#                 "supercont_Wafer8_R2C3_spiral2.3cm_W0.8um_rbw0.1nm_1550-1700nm_vbw100Hz_1mW_2026-02-16-01-56-45.h5"]

# fname_device_hires_Yok = ["Yokogawa_supercont_Wafer8_R2C3_Ring2_W0.8um_gap0.3um_rbw0.1nm_1550-1800nm_6mW_2026-02-15-23-22-21.h5",
#                 "Yokogawa_supercont_Wafer8_R2C3_Ring3_W0.8um_gap0.35um_rbw0.1nm_1550-1800nm_1mW_2026-02-16-00-08-26.h5",
#                 "Yokogawa_supercont_Wafer8_R2C3_Ring4_W0.8um_gap0.4um_rbw0.1nm_1550-1800nm_1mW_2026-02-16-00-28-18.h5",
#                 "Yokogawa_supercont_Wafer8_R2C3_Ring5_W0.8um_gap0.45um_rbw0.1nm_1550-1800nm_1mW_2026-02-16-00-58-35.h5",
#                 ]

# fname_device_coarse = ["supercont_Wafer8_R2C3_Ring0_W0.8um_gap0.2um_600-1700nm_2026-02-15-18-01-35.h5",
#                 "supercont_Wafer8_R2C3_Ring2_W0.8um_gap0.3um_600-1700nm_2026-02-15-18-36-48.h5",
#                 "supercont_Wafer8_R2C3_Ring3_W0.8um_gap0.35um_rbw1nm_600-1700nm_vbw100Hz_6mW_2026-02-15-23-55-14.h5",
#                 "supercont_Wafer8_R2C3_Ring4_W0.8um_gap0.4um_rbw1nm_600-1700nm_vbw100Hz_1mW_2026-02-16-00-39-41.h5",
#                 "supercont_Wafer8_R2C3_Ring5_W0.8um_gap0.45um_rbw1nm_600-1700nm_vbw100Hz_1mW_2026-02-16-00-53-34.h5",
#                 "supercont_Wafer8_R2C3_Ring7_W0.8um_gap0.55um_rbw1nm_600-1700nm_vbw100Hz_1mW_2026-02-16-01-33-09.h5",
#                 "supercont_Wafer8_R2C3_spiral2.3cm_W0.8um_rbw1nm_600-1700nm_vbw100Hz_1mW_2026-02-16-01-58-49.h5"]
# fname_single = ["supercont_Wafer8_R2C3_Ring0_W0.8um_gap0.2um_1100-1400nm_2026-02-15-18-07-55.h5"]
    # "supercont_Wafer8_R2C3_Ring2_W0.8um_gap0.3um_hiRes_1550-1700nm_vbw100Hz_2026-02-15-22-45-59.h5",]
                # "Yokogawa_supercont_Wafer8_R2C3_Ring2_W0.8um_gap0.3um_rbw0.05nm_1550-1800nm_6mW_2026-02-15-23-27-24.h5",]
               
    # "supercont_Wafer8_R2C3_Ring2_W0.8um_gap0.3um_hiRes_1600-1700nm_2026-02-15-18-39-45.h5",
    # "supercont_Wafer8_R2C3_Ring2_W0.8um_gap0.3um_hiRes_1100-1400nm_2026-02-15-18-38-16.h5"]
# dev_label = ["R0_gap0.2um","R2_gap0.3um", "R3_gap0.35um", "R4_gap0.4um", "R5_gap0.45um", "R7_gap0.55um","spiral_2.3cm"]
# dev_label = ["R2_gap0.3um", "R3_gap0.35um", "R4_gap0.4um", "R5_gap0.45um", "R7_gap0.55um","spiral_2.3cm"]

# ----------------------------------------------------------------------------
# 10/21/25 Initial supercontinuum transmission testing of Wafer8 R2C3 Euler Spirals
# sample_dir = 'LL_SiN_EulerSpirals_Test1'
# fname_noise = ["supercont_osa_noiseFloor_2025-10-21-13-37-15.h5"]
# fname_superc = ["supercont_atten_2025-10-17-18-49-10.h5"]
# fname_probe = ["supercont_probe-to-probe_2025-10-21-14-26-53.h5",
#                 "supercont_probe-to-probe_2_2025-10-21-14-34-20.h5",
#                 "supercont_probe-to-probe_3_2025-10-21-14-36-49.h5"]
# fname_device = ["supercont_Wafer8_W0.8um_ref_2025-10-21-13-38-19.h5",
#                 "supercont_Wafer8_W1.4um_ref_2025-10-21-13-51-07.h5",
#                 "supercont_Wafer8_W1.4um_L2.3cm_2025-10-16-23-33-42.h5",
#                 "supercont_Wafer8_R2C3_W2.0um_ref_2025-10-21-13-59-48.h5"]
# dev_label = ["W0.8um_ref","W1.4um_ref", "W1.4um_L2.3cm", "W2.0um_ref"]

# Plot all files in select fname array (initial screening)
fname_arr = fname_probe
fig, ax = plt.subplots(2, 1, figsize=(5,4), tight_layout=True)

x_0, x_f = 600, 1700
for fname in fname_arr: 
    fpath = data_dir / sample_dir / Path(fname)
    ds = load_hdf5(fpath=fpath)
    
    ax[0].plot(ds["wavelength"].to(u.nm), ds["psd"])
    ax[0].set_ylabel('PSD (dBm/nm)')
    ax[1].plot(ds['wavelength'].to(u.nm), ds['amplitude'].to(u.dBm))
    ax[1].set_ylabel('Power (dBm)')
    # ax[1].set_ylim((-90,-60))
    # ax[0].set_ylim((-80,-50))
    ax[1].set_xlabel("Wavelength (nm)")
    # ax[1].set_xlim((x_0, x_f))
    # ax[0].set_xlim((x_0, x_f))
    # plt.legend()

#%% Plot Unormalized PSD >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

ds_noise = load_hdf5(fpath = data_dir/sample_dir/Path(fname_noise[0]))
# ds_superc = load_hdf5(fpath = data_dir/sample_dir/Path(fname_superc[0]))
ds_probe = load_hdf5(fpath = data_dir/sample_dir/Path(fname_probe[0]))

fig, ax = plt.subplots(1, 1, figsize=(5,4), tight_layout=True)
ax.plot(ds_noise["wavelength"].to(u.nm), ds_noise["psd"], color="silver")
# ax.plot(ds_superc["wavelength"].to(u.nm), ds_superc["psd"], color="r", label="source")
# ax.plot(ds_probe["wavelength"].to(u.nm), ds_probe["psd"], color="k", label="probe")
ax.set_ylabel('PSD (dBm/nm)')
ax.set_xlabel("Wavelength (nm)")

x_0, x_f = 600, 1700
for ind, fname in enumerate(fname_device): 
    fpath = data_dir / sample_dir / Path(fname)
    ds = load_hdf5(fpath=fpath)
    
    ax.plot(ds["wavelength"].to(u.nm), ds["psd"], label=dev_label[ind])
    # ax.set_xlim((x_0, x_f))
    # ax.set_ylim((-80, -50))
plt.legend(loc="lower left")

#%% Plot normalized PSD >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#Assuming wavelength axis same for probe-to-probe and device measurments
fig, ax = plt.subplots(1, 1, figsize=(5,4), tight_layout=True)
# ax.plot(ds_noise["wavelength"].to(u.nm), ds_noise["psd"], color="silver")
ymin = -80
ymax = 10
for ind, fname in enumerate(fname_device): 
    fpath = data_dir / sample_dir / Path(fname)
    ds = load_hdf5(fpath=fpath)
    T = ds["psd"] - ds_probe["psd"]
    ax.plot(ds["wavelength"].to(u.nm), T, label=dev_label[ind])
    ax.plot(ds["wavelength"].to(u.nm), ds_noise["psd"] - ds_probe["psd"], color="silver")
    ax.vlines(1000, ymin, ymax, 'k', 'dotted')
    # ax.vlines(1550, ymin, ymax, 'k', 'dotted')
    # ax.vlines(1038, ymin, ymax, 'k', 'dotted')
    # ax.vlines(1080, ymin, ymax, 'k', 'dotted')
    # ax.vlines(1185, ymin, ymax, 'k', 'dotted')
    # ax.vlines(1300, ymin, ymax, 'k', 'dotted')
    ax.set_ylim((ymin,ymax))
ax.set_ylabel('Transmission (dB)')
ax.set_xlabel("Wavelength (nm)")
plt.legend(loc="lower left")
#%% Save Transmission plot >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sname = 'fname'
fig_name = data_dir / Path(fname + '_' + '.png')

fig, ax = plt.subplots(1, 1, figsize=(4,3), tight_layout=True)
ax.plot(wav, T)
ax.set_ylabel('Transmission (dB)')
ax.set_xlabel("Wavelength (nm)")


# %%

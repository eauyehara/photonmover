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

# ------------------------------------------------------------
# 3/11/26 Test BLP01-1550R-25
# sample_dir = "LL_dilute_SiN"
# fname_noise = ["supercont_osa_noiseFloor_2026-03-11-23-44-23.h5"]
# fname_device = [
#                 # "SBS19A_W11_00_10m_600-1700nm_2026-03-11-23-46-28.h5",
#                 # "SBS19A_W11_01_10m_600-1700nm_2026-03-12-00-45-11.h5",
#                 # "SBS19A_W11_10_10m_600-1700nm_pol_2026-03-12-20-32-14.h5",
#                 # "SBS19A_W11_10_10m_600-1700nm_2026-03-12-20-16-55.h5",
#                 # "SBS19A_W8_00_10m_600-1700nm_2026-03-12-22-03-05.h5",
#                 # "SBS19A_W8_01_10m_600-1700nm_2026-03-12-23-07-00.h5",
#                 # "SBS19A_W5_00_1m_600-1700nm_2026-03-13-00-07-55.h5",
#                 # "SBS19A_W5_00_5m_600-1700nm_2026-03-13-01-08-08.h5",
                
#                 # "SBS19A_W5_01_1m_600-1700nm_2026-03-20-20-46-23.h5",
#                 # "SBS19A_W5_01_1m_600-1700nm_DMSO_2026-03-20-21-00-10.h5",
#                 # "SBS19A_W5_01_1m_600-1700nm_DMSO_realign_2026-03-20-21-14-09.h5",
#                 # "SBS19A_W5_01_1m_600-1700nm_DMSO_realign2_2026-03-20-21-33-11.h5",
#                 # "SBS19A_W5_01_5m_600-1700nm_DMSO_2026-03-20-21-53-58.h5",
                
#                 # "SBS19A_W2_00_2m_600-1700nm_TZpolish_2026-03-22-19-28-46.h5",
#                 # "SBS19A_W2_00_3m_600-1700nm_TZpolish_2026-03-22-17-38-49.h5",
#                 # "SBS19A_W2_00_2m_600-1700nm_ImmChamb_dmso_2026-03-25-01-05-50.h5",
#                 # "SBS19A_W2_00_3m_600-1700nm_ImmChamb_dmso_2_2026-03-25-00-46-56.h5",
#                 # # "SBS19A_W2_00_3m_600-1700nm_ImmChamb_dmso_2026-03-25-00-38-21.h5"
#                 # "SBS19A_W2_00_testwg1_600-1700nm_ImmChamb_dmso_2026-03-25-01-00-37.h5",
#                 # "SBS19A_W2_00_testwg2_600-1700nm_ImmChamb_dmso_2026-03-25-00-55-23.h5",
                
#                 "SBS19A_W2_00_testwg2_600-1700nm_2026-03-26-21-09-58.h5",
#                 # "SBS19A_W2_00_testwg2_600-1700nm_1_2026-03-26-21-11-47.h5",
#                 "SBS19A_W2_00_testwg2_600-1700nm_ipa_2026-03-26-21-18-48.h5",
#                 "SBS19A_W2_00_testwg2_600-1700nm_Bcarotene_2026-03-26-21-48-50.h5",
#                 "SBS19A_W2_00_testwg2_600-1700nm_DIwater_2026-03-26-00-15-55.h5",
#                 # "SBS19A_W2_00_3m_600-1700nm_DIwater_2026-03-25-23-59-07.h5",
#                 # "SBS19A_W2_00_3m_600-1700nm_DIwater2_2026-03-26-00-08-00.h5"
#                 ]
# fname_probe = [
#     # "supercont_probe-to-probe_OBJ_2026-03-12-00-57-18.h5",
#     "supercont_probe-to-probe_OBJ_2026-03-12-01-03-41.h5",
#     # "SBS19A_W2_00_testwg2_600-1700nm_TZpolish_2026-03-22-18-32-32.h5",
#     # "SBS19A_W2_00_testwg3_600-1700nm_TZpolish_2026-03-22-19-10-44.h5"
#     ] #["supercont_probe-to-probe_OBJ_2026-03-10-22-53-19.h5"]
# dev_label = ["W2_testwg2", "W2_testwg2_ipa", "W2_testwg2_Bcarot", "W2_testwg2_DI"]
# dev_label = ["W5_1m", "W5_1m_dmso", "W5_1m_dmso1", "W5_1m_dmso2", "W5_1m_dmso3"]
# dev_label = ["W2_2m", "W2_3m", "W2_2m_dmso", "W2_3m_dmso", "W2_testwg1_dmso", "W2_testwg2_dmso"]
# ---------------------------------------------------------------------
# 3/11/26 Tarik Alluxa - objectives
# ------------------------------------------------------------
# 3/11/26 Test BLP01-1550R-25
# sample_dir = "test_BLP01-1550R"
# fname_noise = ["supercont_osa_noiseFloor_2026-03-10-22-51-39.h5"]
# fname_device = ["BLP01_1550R_10deg_600-1700nm_2026-03-11-00-15-36.h5",
#                 "BLP01_1550R_30deg_600-1700nm_2026-03-11-00-12-51.h5",
#                 "BLP01_1550R_45deg_600-1700nm_2026-03-11-00-09-31.h5"]
# fname_probe = ["supercont_probe-to-probe_OBJ_EU_2026-03-11-00-17-09.h5"] #["supercont_probe-to-probe_OBJ_2026-03-10-22-53-19.h5"]
# dev_label = ["10 deg", "30 deg", "45 deg"]

# ---------------------------------------------------------------------
# 3/11/26 Tarik Alluxa - objectives
# sample_dir = "Tarik_Alluxa2"
# fname_noise = ["supercont_osa_noiseFloor_2026-03-10-22-51-39.h5"]
# fname_device = [
#                 # "Alluxa895.8nm_LP830nm_SP100nm_600-1700nm_2026-03-10-22-56-59.h5", # Redundant
#                 "Alluxa895.8nm_LP830nm_SP100nm_600-1700nm_2_2026-03-10-23-08-19.h5",
#                 "LP830nm_SP100nm_600-1700nm_2026-03-10-23-14-39.h5",
#                 # "LP830nm_SP100nm_600-1700nm_2_2026-03-10-23-19-11.h5"  # redundant
#                 # "LP830nm_SP950_SP1000nm_600-1700nm_2026-03-10-23-29-16.h5",  #misaligned
#                 # "LP830nm_SP950_SP1000nm_600-1700nm_2_2026-03-10-23-30-48.h5",  #misalgined
#                 "LP830nm_SP950_600-1700nm_2026-03-10-23-32-44.h5",
#                 # "LP830nm_SP950_600-1700nm_realign_2026-03-10-23-41-57.h5",  # Realign
#                 # "LP830nm_Semrock1000_600-1700nm_realign_2026-03-10-23-47-45.h5"  # Realign
#                 ]
# fname_probe = [
#             "supercont_probe-to-probe_OBJ_2026-03-10-22-53-19.h5",
#             #    "supercont_probe-to-probe_OBJ_2_2026-03-10-23-10-04.h5",
#             #    "supercont_probe-to-probe_OBJ_realign_2026-03-10-23-44-05.h5"  # Realign
#                ] 
# dev_label = ["All filters", "LP830_SP1000", "LP830_SP950"]
# # dev_label = ["SP950", "Semrock1000"]

# ----------------------------------------------------------------------
# 2/23/26 Dilute SiN SBS15A - 0.95 NA objectives
# sample_dir = "LL_dilute_SiN"
# fname_noise = ["supercont_osa_noiseFloor_2026-02-22-23-32-03.h5"]
# fname_device = ["SBS15A_1550_rbw1nm_600-1700nm_vbw400Hz_2026-02-22-23-30-10.h5",
#                 # "SBS15A_1550_AC050-008-C_rbw1nm_600-1700nm_vbw400Hz_2026-02-23-20-59-47.h5",
#                 # "SBS15A_1550_AC050-008-C_FESH0950_rbw1nm_600-1700nm_vbw400Hz_2026-02-23-21-07-43.h5",
#                 # "SBS15A_1550_AC050-008-C_SiPM_rbw1nm_600-1700nm_vbw400Hz_2026-02-23-21-16-51.h5",
#                 # "SBS15A_1550_SiPM_rbw1nm_600-1700nm_vbw400Hz_2026-02-22-23-35-52.h5"
#                 ]
# fname_probe = ["supercont_probe-to-probe_0_2026-02-22-22-10-08.h5"]
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
#     # "supercont_Wafer8_R2C3_Ring2_W0.8um_gap0.3um_hiRes_1550-1700nm_vbw100Hz_2026-02-15-22-45-59.h5",]
#                 # "Yokogawa_supercont_Wafer8_R2C3_Ring2_W0.8um_gap0.3um_rbw0.05nm_1550-1800nm_6mW_2026-02-15-23-27-24.h5",]
                
# fname_device = ["supercont_Wafer8_R2C3_spiral2.3cm_W0.8um_rbw1nm_600-1700nm_vbw100Hz_1mW_2026-02-16-01-58-49.h5"]
               
    # "supercont_Wafer8_R2C3_Ring2_W0.8um_gap0.3um_hiRes_1600-1700nm_2026-02-15-18-39-45.h5",
    # "supercont_Wafer8_R2C3_Ring2_W0.8um_gap0.3um_hiRes_1100-1400nm_2026-02-15-18-38-16.h5"]
# dev_label = ["W8_R2C3_2.3cm"]
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
fname_arr = fname_device
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
    plt.legend()

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
fig, ax0 = plt.subplots(1, 1, figsize=(5,4), tight_layout=True)
# ax0.plot(ds_noise["wavelength"].to(u.nm), ds_noise["psd"], color="silver")
ymin = -70
ymax = -00
for ind, fname in enumerate(fname_device): 
    fpath = data_dir / sample_dir / Path(fname)
    ds = load_hdf5(fpath=fpath)
    T = ds["psd"] - ds_probe["psd"]
    ax0.plot(ds["wavelength"].to(u.nm), T, label=dev_label[ind])
    # ax0.plot(ds["wavelength"].to(u.nm), ds_noise["psd"] - ds_probe["psd"], color="silver")
    # ax0.vlines(1000, ymin, ymax, 'k', 'dotted')
    # ax0.vlines(1550, ymin, ymax, 'k', 'dotted')
    # ax0.vlines(1520, ymin, ymax, 'k', 'dotted')
    # ax0.vlines(1064, ymin, ymax, 'k', 'dotted')
    # ax0.vlines(1383, ymin, ymax, 'k', 'dotted')
    # ax0.vlines(1300, ymin, ymax, 'k', 'dotted')
    ax0.set_ylim((ymin,ymax))
ax0.set_ylabel('Transmission (dB)')
ax0.set_xlabel("Wavelength (nm)")
plt.legend(loc="lower left")
#%% Save Transmission plot >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sname = 'fname'
fig_name = data_dir / Path(fname + '_' + '.png')

fig, ax = plt.subplots(1, 1, figsize=(4,3), tight_layout=True)
ax.plot(wav, T)
ax.set_ylabel('Transmission (dB)')
ax.set_xlabel("Wavelength (nm)")


# %%

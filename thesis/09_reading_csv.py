# -*- coding: utf-8 -*-
# import numpy as np
import pandas as pd
import pathlib as pl
import tikzplotlib as tkz
import matplotlib.pyplot as plt

# directory = "1128_2111_data_77K_p30_r2"
# directory = "1206_0720_data_77K_p30_r2_plots_second_guess"
directory = "v2"


#%% --- voltages ---

filename = "HB_Ui"
path = pl.Path(__file__).parent / "visualize" / directory / (filename + ".csv")

columns = ["Time", "U", "V", "W"]
df = pd.read_csv(path, usecols=columns) # df = dataframe



fig = plt.figure(dpi=300, figsize=(10,7))
ax1 = plt.subplot()
plt.grid()

ax1.set_xlabel("Time")
ax1.set_ylabel("Voltage")
ax1.set_xlim(4.5018, 5.40216)
ax1.set_ylim(-690, 690)
ax1.plot(df.Time, df.U, color = "blue", label="U")
ax1.plot(df.Time, df.V, color = "red", label="V")
ax1.plot(df.Time, df.W, color = "green", label="W")
# ax1.plot([1.4406, 1.6807], [563.38]*2)
# ax1.plot([1.4406, 1.6807], [-563.38]*2)
ax1.legend()

tkz.clean_figure()
# tkz.save("aus_comsol_Ui_HB_30K_p8.tex")
tkz.save("aus_comsol_Ui_HB_30K_p8_second_guess.tex")

# tkz.save("aus_comsol_Ui_HB_77K_p30.tex")
# tkz.save("aus_comsol_Ui_HB_77K_p30_second_guess.tex")


#%% --- current densities ---

filename = "HB_J"
path = pl.Path(__file__).parent / "visualize" / directory / (filename + ".csv")

columns = ["Time", "U", "V", "W"]
df = pd.read_csv(path, usecols=columns) # df = dataframe

fig = plt.figure(dpi=300, figsize=(10,7))
ax1 = plt.subplot()
plt.grid()

ax1.set_xlabel("Time")
ax1.set_ylabel("J")
ax1.set_xlim(4.5018, 5.40216) # ax1.set_xlim(1.4406, 1.6807)
ax1.set_ylim(-465, 465) # ax1.set_ylim(-60, 60)
fac = 1e-6
ax1.plot(df.Time, df.U*fac, color = "blue", label="U")
ax1.plot(df.Time, df.V*fac, color = "red", label="V")
ax1.plot(df.Time, df.W*fac, color = "green", label="W")
ax1.legend()

tkz.clean_figure()
# tkz.save("aus_comsol_Je_HB_30K_p8.tex")
tkz.save("aus_comsol_Je_HB_30K_p8_second_guess.tex")

# tkz.save("aus_comsol_Je_HB_77K_p30.tex")
# tkz.save("aus_comsol_Je_HB_77K_p30_second_guess.tex")

#%% --- torque ---

filename = "HB_M"
path = pl.Path(__file__).parent / "visualize" / directory / (filename + ".csv")

columns = ["Time", "M"]
df = pd.read_csv(path, usecols=columns) # df = dataframe

fig = plt.figure(dpi=300, figsize=(10,7))
ax1 = plt.subplot()
plt.grid()

ax1.set_xlabel("Time")
ax1.set_ylabel("M")
ax1.set_xlim(4.5018, 5.40216) # ax1.set_xlim(1.4406, 1.6807)
# ax1.set_ylim(-690, 690)
ax1.set_ylim(6.08, 6.16) # ax1.set_ylim(3.52, 3.6) # ax1.set_ylim(4.3, 4.35)
ax1.plot(df.Time, [M*1e-6 for M in df.M], color = "blue", label="M")

tkz.clean_figure()
# tkz.save("aus_comsol_Me_HB_30K_p8.tex")
tkz.save("aus_comsol_Me_HB_30K_p8_second_guess.tex")

# tkz.save("aus_comsol_Me_HB_77K_p30.tex")
# tkz.save("aus_comsol_Me_HB_77K_p30_second_guess.tex")

# filename = "HB_torque_via_maxwell"
# path = pl.Path(__file__).parent / "visualize" / directory / (filename + ".csv")

# columns = ["Time", "M"]
# df = pd.read_csv(path, usecols=columns) # df = dataframe
# ax1.plot(df.Time, df.M, color = "red", label="M")

#%% -- isovac data ---

directory = "isovac"

filename = "isovac_BH"
path = pl.Path(__file__).parent / "visualize" / directory / (filename + ".csv")

columns = ["B", "H"]
df = pd.read_csv(path, usecols=columns) # df = dataframe

fig = plt.figure(dpi=300, figsize=(10,7))
ax1 = plt.subplot()
plt.grid()

ax1.set_xlabel("B")
ax1.set_ylabel("H")
ax1.set_xlim(0,2)
ax1.set_ylim(0, 500)
ax1.plot(df.B, [M for M in df.H], color = "blue", label="M")

#%% -- isovac data ---
import matplotlib.ticker as mticker

directory = "isovac"

filename = "isovac_HB"
path = pl.Path(__file__).parent / "visualize" / directory / (filename + ".csv")

columns = ["H", "B"]
df = pd.read_csv(path, usecols=columns) # df = dataframe

fig = plt.figure(dpi=300, figsize=(7,10))
ax1 = plt.subplot()
plt.grid(True, which="both")

ax1.set_xlabel("H")
ax1.set_ylabel("B")
ax1.set_xscale('log')
ax1.xaxis.set_major_formatter(mticker.ScalarFormatter())
ax1.set_xlim(10,10000)
ax1.set_ylim(0, 1.8)
ax1.plot(df.H, [M for M in df.B], color = "blue", label="M")

tkz.clean_figure()
tkz.save("app_isovac_HB_edit.tex")

#%% -- isovac data ---
import matplotlib.ticker as mticker

directory = "isovac"

filename = "isovac_HB"
path = pl.Path(__file__).parent / "visualize" / directory / (filename + ".csv")

columns = ["H", "B"]
df = pd.read_csv(path, usecols=columns) # df = dataframe

fig = plt.figure(dpi=300, figsize=(7,10))
ax1 = plt.subplot()
plt.grid(True, which="both")

ax1.set_xlabel("H")
ax1.set_ylabel("B")
ax1.set_xscale('log')
ax1.xaxis.set_major_formatter(mticker.ScalarFormatter())
ax1.set_xlim(10,2000000)
ax1.set_ylim(0, 4)
ax1.plot(df.H, [M for M in df.B], color = "blue", label="M")

tkz.clean_figure()
# tkz.save("app_isovac_HB_edit.tex")


#%% -- isovac data ---
import matplotlib.ticker as mticker

directory = "isovac"

filename = "isovac_BH"
path = pl.Path(__file__).parent / "visualize" / directory / (filename + ".csv")

columns = ["B", "H"]
df = pd.read_csv(path, usecols=columns) # df = dataframe

fig = plt.figure(dpi=300, figsize=(7,10))
ax1 = plt.subplot()
plt.grid(True, which="both")

ax1.set_xlabel("B")
ax1.set_ylabel("H")
ax1.set_yscale('log')
# ax1.xaxis.set_major_formatter(mticker.ScalarFormatter())
ax1.set_xlim(0,6)
ax1.set_ylim(10, 3e6)
ax1.plot(df.B, [M for M in df.H], color = "blue", label="M")

tkz.clean_figure()
tkz.save("app_isovac_BH_edit.tex")
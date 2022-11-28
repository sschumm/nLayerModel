# -*- coding: utf-8 -*-
# import numpy as np
import pandas as pd
import pathlib as pl
import tikzplotlib as tkz
import matplotlib.pyplot as plt

directory = "1128_2111_data_77K_p30_r2"


#%% --- voltages ---

filename = "HB_induced_voltages"
path = pl.Path(__file__).parent / "visualize" / directory / (filename + ".csv")

columns = ["Time", "U", "V", "W"]
df = pd.read_csv(path, usecols=columns) # df = dataframe



fig = plt.figure(dpi=300, figsize=(10,7))
ax1 = plt.subplot()
plt.grid()

ax1.set_xlabel("Time")
ax1.set_ylabel("Voltage")
ax1.plot(df.Time, df.U, color = "blue", label="U")
ax1.plot(df.Time, df.V, color = "red", label="V")
ax1.plot(df.Time, df.W, color = "green", label="W")
# ax1.legend()

tkz.clean_figure()
# tkz.save("aus_Ui_HB_77K_p30.tex")


#%% --- current densities ---

filename = "HB_current_densities"
path = pl.Path(__file__).parent / "visualize" / directory / (filename + ".csv")

columns = ["Time", "U", "V", "W"]
df = pd.read_csv(path, usecols=columns) # df = dataframe

fig = plt.figure(dpi=300, figsize=(10,7))
ax1 = plt.subplot()
plt.grid()

ax1.set_xlabel("Time")
ax1.set_ylabel("J")
ax1.plot(df.Time, df.U, color = "blue", label="U")
ax1.plot(df.Time, df.V, color = "red", label="V")
ax1.plot(df.Time, df.W, color = "green", label="W")
# ax1.legend()

#%% --- torque ---

filename = "HB_torque_force_calculation"
path = pl.Path(__file__).parent / "visualize" / directory / (filename + ".csv")

columns = ["Time", "M"]
df = pd.read_csv(path, usecols=columns) # df = dataframe

fig = plt.figure(dpi=300, figsize=(10,7))
ax1 = plt.subplot()
plt.grid()

ax1.set_xlabel("Time")
ax1.set_ylabel("M")
ax1.plot(df.Time, df.M, color = "blue", label="M")

filename = "HB_torque_via_maxwell"
path = pl.Path(__file__).parent / "visualize" / directory / (filename + ".csv")

columns = ["Time", "M"]
df = pd.read_csv(path, usecols=columns) # df = dataframe
ax1.plot(df.Time, df.M, color = "red", label="M")

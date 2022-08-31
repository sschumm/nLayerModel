# -*- coding: utf-8 -*-
import numpy as np
np.set_printoptions(suppress=True, linewidth=250, precision=2)

from scipy.constants import pi
from modules import Model, MagneticLayer, AirLayer, CurrentLoading
from modules.plot import RadialMultiPlot, PlanePlot
from analytics.precalcs import kb, kd, kp, K, taup
from data import Generator as gn
from data import StatorWinding_Cu as sw
from data import FieldWinding as fw


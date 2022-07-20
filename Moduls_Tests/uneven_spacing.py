# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 09:49:26 2022

@author: svens
"""

# from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
import Moduls_Tests.somemath as sm
from scipy.interpolate import griddata

xx = np.array([[-32.77352506, -32.50517324, -32.30341846, -32.12060867, -31.99103968],
   [-32.88670112, -32.63078693, -32.42892793, -32.2527705 , -32.11911059],
   [-32.99884749, -32.75419286, -32.55377179, -32.38220417,-32.24664094],
   [-33.10888993, -32.87495033, -32.67707405, -32.50885765,-32.37311971],
   [-33.2179889 , -32.99317728, -32.79848554, -32.63304009,-32.49815164],
   [-33.32651265, -33.10917094, -32.91791567, -32.75496914,-32.62146004],
   [-33.43449219, -33.22321261, -33.03537769, -32.87477271,-32.74286757],
   [-33.54184645, -33.33551502, -33.15092328, -32.99252837,-32.86227065],
   [-33.64847256, -33.44622558, -33.26461466, -33.1082882 ,-32.97961644],
   [-33.7542787 , -33.55544297, -33.37651245, -33.22209182,-33.09488502],
   [-33.85919362, -33.66323316, -33.48667092, -33.33397251,-33.20807667],
   [-33.96316738, -33.76964127, -33.59513666, -33.44395994,-33.31920308],
   [-34.06616899, -33.87469962, -33.70194888, -33.55208114,-33.42828168],
   [-34.16818354, -33.97843289, -33.80714029, -33.6583608 , -33.5353318 ],
   [-34.26920936, -34.08086103, -33.91073804, -33.76282129,-33.64037232],
   [-34.36925566, -34.18200076, -34.01276449, -33.86548266,-33.74342012],
   [-34.46834041, -34.28186594, -34.11323791, -33.96636282,-33.84448915],
   [-34.56648851, -34.38046701, -34.21217281, -34.06547784,-33.94358994],
   [-34.66372971, -34.4778096 , -34.30958029, -34.16284263, -34.0407292 ],
   [-34.76009619, -34.57389229, -34.40546833, -34.2584719 ,-34.13590974],
   [-34.85561889, -34.66870332, -34.49984222, -34.35238178,-34.22913045],
   [-34.95032165, -34.76221629, -34.59270529, -34.44459209,-34.32038648],
   [-35.04421102, -34.85438466, -34.6840602 , -34.5351297 ,-34.40966975],
   [-35.13725786, -34.94513595, -34.77391107, -34.62403315,-34.49697031],
   [-35.22936336, -35.03436723, -34.86226647, -34.7113589 , -34.5822793 ],
   [-35.3202956 , -35.12194712, -34.9491432 , -34.79718943,-34.66559567],
   [-35.40957149, -35.20773648, -35.03457067, -34.88164276,-34.74693998],
   [-35.49624617, -35.29165455, -35.11859869, -34.96488179,-34.82637986],
   [-35.57858891, -35.3738341 , -35.20134195, -35.04711944,-34.90406862],
   [-35.65407127, -35.45484625, -35.28327484, -35.12861303,-34.98028311]])


yy = np.array([[-11.3529916 , -10.83017948, -10.36062676,  -9.85499224,-9.36742115],
   [-11.24914312, -10.77486528, -10.30767657,  -9.81790781,-9.33347811],
   [-11.16896123, -10.71827884, -10.25654788,  -9.77873607,-9.29985941],
   [-11.09864806, -10.66157581, -10.20688907,  -9.7389486 ,-9.26669158],
   [-11.03379175, -10.60554773, -10.15836065,  -9.69919353,  -9.2340536 ],
   [-10.97234269, -10.55056628, -10.11076275,  -9.65973621,-9.20197376],
   [-10.91318741, -10.49672964, -10.06398502,  -9.62068888,-9.17044015],
   [-10.85567682, -10.44399733, -10.01795592,  -9.58209354,-9.13941538],
   [-10.79941292, -10.39227109,  -9.97261586,  -9.54395471,-9.10884902],
   [-10.74413933, -10.3414354 ,  -9.92790608,  -9.50625469,-9.07868586],
   [-10.68968138, -10.29137538,  -9.88376528,  -9.46896173,-9.04887037],
   [-10.63591212, -10.24198327,  -9.84012944,  -9.43203502,-9.01934878],
   [-10.58273238, -10.19315956,  -9.79693259,  -9.39542787,-8.99006955],
   [-10.53005851, -10.14481189,  -9.75410767,  -9.35908973,-8.96098324],
   [-10.47781454, -10.09685326,  -9.71158727,  -9.32296761,-8.93204187],
   [-10.42592671, -10.04920008,  -9.66930413,  -9.28700687,-8.90319822],
   [-10.37431922, -10.00177041,  -9.62719161,  -9.25115173,  -8.874405  ],
   [-10.32291042,  -9.95448254,  -9.5851842 ,  -9.21534553,-8.84561388],
   [-10.27160885,  -9.90725399,  -9.54321824,  -9.17953081,  -8.8167744 ],
   [-10.22030865,  -9.86000117,  -9.50123304,  -9.1436493 ,-8.78783251],
   [-10.16888358,  -9.81263985,  -9.45917264,  -9.10764179,-8.75872859],
   [-10.11717908,  -9.7650871 ,  -9.41698825,  -9.07144794,-8.72939466],
   [-10.06500151,  -9.71726547,  -9.37464171,  -9.03500603,-8.69975035],
   [-10.01210356,  -9.66911062,  -9.33211024,  -8.99825263,-8.66969694],
   [ -9.95816574,  -9.62058406,  -9.28939255,  -8.96112225,  -8.639109  ],
   [ -9.90277604,  -9.57169217,  -9.24651752,  -8.92354679,-8.60782312],
   [ -9.84541584,  -9.52250952,  -9.2035585 ,  -8.88545441,-8.57562469],
   [ -9.78546516,  -9.4731919 ,  -9.1606645 ,  -8.84676681,-8.54223813],
   [ -9.72217224,  -9.42392909,  -9.11814094,  -8.80739381,-8.50733329],
   [ -9.65407127,  -9.37474785,  -9.07654968,  -8.76722606,-8.47056622]])

np.set_printoptions(suppress=True, linewidth=200, precision=5)

p = 2
aj = 4
bj = 3

r_i, r_o = 0.05, 3
r = np.linspace(r_i, r_o, 50)
t = np.linspace(0, 2*np.pi, 30)
R, T = np.meshgrid(r, t)

Br = sm.Br_no_k(p, R, T, aj, bj)
Bt = sm.B0_no_k(p, R, T, aj, bj)

X, Y = sm.r0_to_xy(R, T)
U, V = sm.BrB0_to_UV(Br, Bt, T)
speed = np.sqrt(U**2, V**2)

# def streams(ax,xx,yy,u,v,base_map=False):
    
x = np.linspace(X.min(), X.max(), 500)
y = np.linspace(Y.min(), Y.max(), 500)
 
xi, yi = np.meshgrid(x,y)

#then, interpolate your data onto this grid:

px = X.flatten()
py = Y.flatten()
pu = U.flatten()
pv = V.flatten()
pspeed = speed.flatten()

# points = np.stack((px, py), axis = -1)
# xi = np.stack((x, y), axis = -1)
# gu = griddata(points=points, values=pu, xi=xi)

gu = griddata((px,py), pu, (xi,yi))
gv = griddata((px,py), pv, (xi,yi))
gspeed = griddata((px,py), pspeed, (xi,yi))

# =============================================================================
# gu = griddata(zip(px,py), pu, (xi,yi))
# gv = griddata(zip(px,py), pv, (xi,yi))
# gspeed = griddata(zip(px,py), pspeed, (xi,yi))
# =============================================================================


plt.figure(figsize=(10, 10))
ax = plt.subplot()

# ax.contour(xx,yy,speed, colors='k', alpha=0.4)
# ax.plot(xx,yy,'-k',alpha=0.3)
# ax.plot(xx.T,yy.T,'-k',alpha=0.3)
# ax.plot(xi,yi,'-b',alpha=0.1)
# ax.plot(xi.T,yi.T,'-b',alpha=0.1)

# b = ax.quiver(X, Y, U, V)
c = ax.streamplot(x,y,gu,gv, density=2, 
                  linewidth=1, 
                  color=gspeed, 
                  cmap=plt.cm.jet)

# =============================================================================
# u = 1/2*np.sin(xx)*-3*np.cos(yy**2)
# v = 2*np.sin(xx)*3*np.cos(yy)
# speed = np.sqrt((u**2)+(v**2))
# 
# plt.ion()
# fig,(ax1,ax2) = plt.subplots(1,2)
# 
# m = Basemap(llcrnrlon=xx.min(),llcrnrlat=yy.min(),
# urcrnrlon=xx.max(), urcrnrlat=yy.max(),
# projection='merc',resolution='i',ax=ax1)
# 
# m.drawcoastlines()
# m.drawcountries()
# m.fillcontinents(color='0.86')  
# m.contourf(xx,yy,speed,latlon=True)
# m.plot(xx,yy,'-k',alpha=0.3,latlon=True)
# m.plot(xx.T,yy.T,'-k',alpha=0.3,latlon=True)
# m.quiver(xx,yy,u,v,latlon=True)
# 
# 
# m2 = Basemap(llcrnrlon=xx.min(),llcrnrlat=yy.min(),
# urcrnrlon=xx.max(), urcrnrlat=yy.max(),
# projection='merc',resolution='i',ax=ax2)
# 
# m2.drawcoastlines()
# m2.drawcountries()
# m2.fillcontinents(color='0.86') 
# m2.contourf(xx,yy,speed,latlon=True)
# m2.plot(xx,yy,'-k',alpha=0.3,latlon=True)
# m2.plot(xx.T,yy.T,'-k',alpha=0.3,latlon=True)
# m2.streamplot(xx,yy,u,v)
# =============================================================================


#%%
import numpy as np
import matplotlib.pyplot as plt
import Moduls_Tests.somemath as sm
from scipy.interpolate import griddata

x = np.linspace(1, 7, 10)
y = np.array([5, 2, 3, 1., 9, 2, 4, 1, 2, 3])
X, Y = np.meshgrid(x, y)

U = np.sin(Y) * X
V = X * (Y - 3)

#####################################

x = np.linspace(X.min(), X.max(), 5)
y = np.linspace(Y.min(), Y.max(), 5)

xi, yi = np.meshgrid(x,y)

px = X.flatten()
py = Y.flatten()
pu = U.flatten()
pv = V.flatten()
# pspeed = speed.flatten()

gu = griddata((px,py), pu, (xi,yi))

# plt.figure(figsize=(10, 10))
# ax = plt.subplot()

# b = ax.quiver(X, Y, U, V)
# c = ax.streamplot(X,Y,U,V, density=2, linewidth=1, color="b", cmap=plt.cm.jet)




















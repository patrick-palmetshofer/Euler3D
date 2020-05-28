import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

deg = "20"
it = ""
#1041.566000 555.500000 867.970000
vel = "555.500000"

limiter = "with"

gridpath = "mesh/Grid"+deg+"deg.grd"
#solpath = "raw/"+limiter+"limiter/Grid"+deg+"deg"+vel+".res"

fname = "solution/test1000"

solpath = fname+".sol"

with open(gridpath) as f:
    first_line = f.readline()
	
print(first_line)

with open(gridpath) as f:
    first_line = f.readline()
	
print(first_line)

linevec = first_line.split()

nxi = int(linevec[0])
neta = int(linevec[1])
	


ncells = [int(vn) for vn in first_line.split()]

grid = pd.read_csv(gridpath,skiprows=1,header=None,delimiter="\s")
solution = pd.read_csv(solpath,skiprows=1,delimiter=",\t")



X =  grid.loc[:,0].values.reshape((nxi,neta))
Y =  grid.loc[:,1].values.reshape((nxi,neta))
data=[0]*4
for i in range(4):
    data[i] = solution.iloc[:,i].values.reshape((nxi-1,neta-1))

fig1, ax = plt.subplots(1,1)
#fig1.set_figheight(15)
#fig1.set_figwidth(22)

cmap = "plasma"
R = 1005*(1-1/1.4)
gamma = 1.4
Ma = np.sqrt((data[1]**2+data[2]**2)/(gamma*R*data[3]))
Ma0 = Ma[0,0]
d = Ma
plot = ax.pcolormesh(X,Y,d,cmap=cmap,edgecolors="black")
ax.axis('equal')
fig1.suptitle("phi="+str(deg)+"°, Ma="+str(round(Ma0,1)))
#fig1.g("GridMach"+limiter+"Limiter"+deg+"deg"+vel+"vel"+".pgf")


fig1, axes = plt.subplots(2,2)
axes = axes.reshape((4,1))
#fig1.set_figheight(15)
#fig1.set_figwidth(22)

cmap = "plasma"
titles = ["Density [kg/m^3]", "u [m/s]", "v [m/s]", "T [K]"]
R = 1005*(1-1/1.4)
gamma = 1.4
Ma = np.sqrt((data[1]**2+data[2]**2)/(gamma*R*data[3]))
Ma0 = Ma[0,0]
for i,ax in enumerate(axes):
    d = data[i]
    ax = axes[i,0]
    plot = ax.pcolormesh(X,Y,d,cmap=cmap)
    ax.axis('equal')
    fig1.colorbar(plot, ax = ax)
    ax.set_title(titles[i])
fig1.suptitle("phi="+str(deg)+"°, Ma="+str(round(Ma0,1)))

# ADD CORRECT CELL MIDPOINTS
Xa = (X[1:,:] + X[:-1,:]) / 2
Xcent = (Xa[:,1:] + Xa[:,:-1]) / 2
Ya = (Y[1:,:] + Y[:-1,:]) / 2
Ycent = (Ya[:,1:] + Ya[:,:-1]) / 2
fig2, axes2 = plt.subplots(2,2)
axes2 = axes2.reshape((4,1))
#fig2.set_figheight(15)
#fig2.set_figwidth(22)
for i,ax in enumerate(axes2):
    d = data[i]
    ax = axes2[i,0]
    plot = ax.contourf(Xcent,Ycent,d,cmap=cmap)
    ax.axis('equal')
    fig2.colorbar(plot, ax = ax)
    ax.set_title(titles[i])
fig2.suptitle("phi="+str(deg)+"°, Ma="+str(round(Ma0,1)))

R = 1005*(1-1/1.4)
gamma = 1.4
Ma = np.sqrt((data[1]**2+data[2]**2)/(gamma*R*data[3]))
figMa, axMa = plt.subplots()
           # figMa.set_figheight(15)
           # figMa.set_figwidth(22)
axMa.axis('equal')
plot = axMa.contourf(Xcent,Ycent,Ma,cmap=cmap)
figMa.colorbar(plot, ax = axMa)
axMa.set_title("Mach number [-] phi="+str(deg)+"°, Ma="+str(round(Ma0,1)))
#figMa.savefig("Ma"+limiter+"Limiter"+deg+"deg"+vel+"vel"+".pgf")

pt = data[0]*(R*data[3]+0.5*(data[1]**2+data[2]**2))/1.01325e5
figpt, axpt = plt.subplots()
           # figMa.set_figheight(15)
           # figMa.set_figwidth(22)
axpt.axis('equal')
plot = axpt.contourf(Xcent,Ycent,pt,cmap=cmap)
figpt.colorbar(plot, ax = axpt)
axpt.set_title("Total pressure [atm] phi="+str(deg)+"°, Ma="+str(round(Ma0,1)))
#figpt.savefig("TotalPressure"+limiter+"Limiter"+deg+"deg"+vel+"vel"+".pgf")

Tt = data[3]+0.5*data[0]/1005*(data[1]**2+data[2]**2)
figTt, axTt = plt.subplots()
           # figMa.set_figheight(15)
           # figMa.set_figwidth(22)
axTt.axis('equal')
plot = axTt.contourf(Xcent,Ycent,Tt,cmap=cmap)
figTt.colorbar(plot, ax = axTt)
axTt.set_title("Total Temperature [K] phi="+str(deg)+"°, Ma="+str(round(Ma0,1)))
#figTt.savefig("TotalTemp"+limiter+"Limiter"+deg+"deg"+vel+"vel"+".pgf")

residualpath = fname+"_Linfty.res"
residuals = pd.read_csv(residualpath,delimiter=",\t")
figRes, axRes = plt.subplots()
#figRes.set_figheight(15)
#figRes.set_figwidth(22)
for col in residuals.columns:
    axRes.plot(residuals.index,residuals[col],label=col)
axRes.set_yscale("log")
#axRes.set_title("Normalized Residual $L_\infty$  phi="+str(deg)+"°, Ma="+str(round(Ma0,1)))
axRes.set_title("Normalized Residual $L_\infty$ ")
axRes.set_xlabel("$L_\infty$")
axRes.set_xlabel("Iteration")
axRes.legend()
#axRes.axis([0,10000,1e-12,1])
axRes.grid()
#figRes.savefig("ResidualsLInf"+limiter+"Limiter"+deg+"deg"+vel+"vel"+".pgf")

residualpath = fname+"_L2.res"
residuals = pd.read_csv(residualpath,delimiter=",\t")
figRes, axRes = plt.subplots()
#figRes.set_figheight(15)
#figRes.set_figwidth(22)
for col in residuals.columns:
    axRes.plot(residuals.index,residuals[col],label=col)
axRes.set_yscale("log")
#axRes.set_title("Normalized Residual $L_2$  phi="+str(deg)+"°, Ma="+str(round(Ma0,1)))
axRes.set_title("Normalized Residual $L_2$")
axRes.set_xlabel("$L_2$")
axRes.set_xlabel("Iteration")
axRes.legend()
#axRes.axis([0,10000,1e-12,1])
axRes.grid()
#figRes.savefig("ResidualsL2"+limiter+"Limiter"+deg+"deg"+vel+"vel"+".pgf")
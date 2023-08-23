import numpy as np
import matplotlib.pyplot as plt
import math
import os
from xlwt import Workbook
from scipy.optimize import curve_fit

CfcalKoef = 0.913 # linear koeficient of proton drop point calibration curve

def read(loadName):
    highChannel, lowChannel, highy, lowy, highGain, lowGain, le, gain = [],[],[],[],[],[],[],[]
    with open(loadName, "r") as f:
        data = f.readlines()
        start,k,h,l=0,0,0,0
        for line in range(len(data)):
            if data[line].find("Counts")>-1:
                start=1
            if start==1:
                words = data[line].split()
                if len(words)==11:
                    highChannel.append(int(words[0]))
                    lowChannel.append(int(words[0]))
                    highy.append(int(words[0])*0.1)
                    lowy.append(int(words[0])*1.5)
                    highGain.append(int(words[10]))
                    lowGain.append(int(words[4]))
                if len(words)==5:
                    lowChannel.append(int(words[0]))
                    lowy.append(int(words[0])*1.5)
                    lowGain.append(int(words[4]))
        for i in range(0, 1262):
            if i<=255:
                le.append(highy[i])
                gain.append(highGain[i])
            if i>255:
                le.append(lowy[i-238])
                gain.append(lowGain[i-238])
    return list(highChannel),list(lowChannel),list(highy),list(lowy), list(highGain),list(lowGain),list(le),list(gain)

wb = Workbook()
normFrequencyConstant,normDoseConstant,center=[],[],[]
name=input("Filename input: ")
fileName=name+".SP2"
xlsName=name+".xls"

highChannel, lowChannel, highy, lowy, highGain, lowGain, linBins, linCount = read(fileName)

def bin(l1, l2, log1, log2, a, b):
    k,Slog,S = 0, 0, 1
    if (a == 0 and b == 0):
        Slog,S = 0, 1
    elif (l1 < log2 and l1 >= log1 and l2 >= log2):
        delta = log2 - l1
        linDelta = l2 - l1
        Slog = b * delta + (a * (delta**2)) / 2
        S = b * linDelta + (a * (linDelta**2)) / 2
    elif l1 < log1 and l2 > log1 and l2 <= log2:
        delta = l2 - log1
        linDelta = l2 - l1
        Slog = b * delta + (a * (delta**2)) / 2
        S = b * linDelta + (a * (linDelta**2)) / 2
    elif (l1 < log1 and l2 > log2):
        delta = log2 - log1
        linDelta = l2 - l1
        Slog = b * delta + (a * (delta**2)) / 2
        S = b * linDelta + (a * (linDelta**2)) / 2
    elif l1 > log1 and l2 < log2:
        Slog, S = 1, 1
    k = Slog / S
    return k

logBins=np.logspace(np.log10(0.1),np.log10(1000),num=51)
for i in range(len(logBins)-1):
    center.append(math.pow(10,np.log10(logBins[i])+(np.log10(logBins[i+1])-np.log10(logBins[i]))/2))
widths=(logBins[1:]-logBins[:-1])
logCount=[[0 for j in range(len(logBins)-1)] for i in range(3)]
event=[[0 for j in range(len(logBins)-1)] for i in range(3)]
# Binning
for i in range(len(logBins)-1):
    for j in range(len(linBins)-1):
        if linBins[j+1]<logBins[i]:
            continue
        if linBins[j]>logBins[i+1]:
            break
        if linBins[j+1]>logBins[i] and linBins[j+1]<=logBins[i+1]:
            logCount[0][i]=logCount[0][i]+linCount[j]
            event[0][i]+=1
        a=linCount[j+1]-linCount[j]
        b=linCount[j]
        q1=bin(linBins[j],linBins[j+1],logBins[i],logBins[i+1],0,b)
        q2=bin(linBins[j],linBins[j+1],logBins[i],logBins[i+1],a,b)
        if q1>0:
            logCount[1][i]=logCount[1][i]+linCount[j]*q1
            event[1][i]+=1
        if q2>0:
            logCount[2][i]=logCount[2][i]+linCount[j]*q2
            event[2][i]+=1
# Normovani a kalkulovani
frequencyDistribution=[[0 for j in range(len(logBins)-1)] for i in range(3)]
doseDistribution=[[0 for j in range(len(logBins)-1)] for i in range(3)]
Ny_Wy=[[0 for j in range(len(logBins)-1)] for i in range(3)]
Ny_Wy_y=[[0 for j in range(len(logBins)-1)] for i in range(3)]
fy=[[0 for j in range(len(logBins)-1)] for i in range(3)]
y_fy=[[0 for j in range(len(logBins)-1)] for i in range(3)]
y_fy_Wy=[[0 for j in range(len(logBins)-1)] for i in range(3)]
dy=[[0 for j in range(len(logBins)-1)] for i in range(3)]
y_dy=[[0 for j in range(len(logBins)-1)] for i in range(3)]
y_dy_Wy=[[0 for j in range(len(logBins)-1)] for i in range(3)]
for j in range(3):
    for i in range(len(logBins)-1):
        if event[j][i]==0:
            event[j][i]=1
        logCount[j][i]=logCount[j][i]/event[j][i]
        Ny_Wy[j][i]=logCount[j][i]*widths[i]
        Ny_Wy_y[j][i]=logCount[j][i]*widths[i]*center[i]

for j in range(3):
    for i in range(len(logBins)-1):
        fy[j][i]=logCount[j][i]/float(sum(Ny_Wy[j]))
        y_fy[j][i]=fy[j][i]*center[i]
        y_fy_Wy[j][i]=center[i]*fy[j][i]*widths[i]
        dy[j][i]=logCount[j][i]*center[i]/float(sum(Ny_Wy_y[j]))
        y_dy[j][i]=dy[j][i]*center[i]
        y_dy_Wy[j][i]=dy[j][i]*center[i]*widths[i]

# select region of p-drop point
def select_p_drop_point(lowChannel,lowGain):
    CfCalChn, CfCalCoun = [],[]
    for i in range(len(lowChannel)):
        if (lowChannel[i] >= 50) and (lowChannel[i] <= 200):
            CfCalChn.append(lowChannel[i])
            CfCalCoun.append(lowGain[i]*lowChannel[i]*lowChannel[i])
    return list(CfCalChn),list(CfCalCoun)

def Fermi(x, a, b, c):
    return a/(1+np.exp(b*(x-c)))

# find first aproximation of p-edge
def find_aprox(CfCalChn,CfCalCoun):
    def Fermi(x, a, b, c):
        return a/(1+np.exp(b*(x-c)))
    p0 = [100000, 0.08, 100]
    CfCalPar = curve_fit(Fermi, CfCalChn, CfCalCoun, p0[:])
    [a, b, c] = CfCalPar[0]
    CfCalChn2, CfCalCoun2 = [],[]
    aaa = int(c-3/abs(b))
    bbb = int(c+5/abs(b))
    for i in range(len(CfCalChn)):
        if (CfCalChn[i] >= aaa) and (CfCalChn[i] <= bbb):
            CfCalChn2.append(CfCalChn[i])
            CfCalCoun2.append(CfCalCoun[i])
#    def background(x, d):
#        return d+0*x
#    param_bounds=([-np.inf,bbb,-np.inf],[np.inf,1024,np.inf])
#    backgroundPAR = curve_fit(Fermi, CfCalChn, CfCalCoun, bounds=param_bounds)
#    d = backgroundPAR[0]
#    print(backgroundPAR[0])
    return list(CfCalChn2),list(CfCalCoun2)

# find parameters of Fermi fit
def find_Fermi(CfCalChn,CfCalCoun):
    def Fermi(x, a, b, c):
        return a/(1+np.exp(b*(x-c)))
    p0 = [100000, 0.08, 100]
    CfCalPar = curve_fit(Fermi, CfCalChn, CfCalCoun, p0[:])
    [a, b, c] = CfCalPar[0]
    return [a, b, c]

# propose approximate change of voltage
def propose_voltage_change(CfcalKoef,c):
    CfpropVolChange = 0
    if abs(c-100) > 3:
        CfpropVolChange = (c-100)/(CfcalKoef)
    return CfpropVolChange

CfCalChn, CfCalCoun = select_p_drop_point(lowChannel,lowGain)
CfCalChn2, CfCalCoun2 = find_aprox(CfCalChn,CfCalCoun)
[a, b, c] = find_Fermi(CfCalChn2,CfCalCoun2)
CfpropVolChange = propose_voltage_change(CfcalKoef,c)

# Ulozeni do excelovskych souboru
sheet1=wb.add_sheet("Counts")
sheet2=wb.add_sheet("Frequency distribution")
sheet3=wb.add_sheet("Dose distribution")
sheet4=wb.add_sheet("Equivalent dose distribution")
sheet5=wb.add_sheet("Calibrations")

sheet1.write(0,0,"High Gain Spectrum")
sheet1.write(1,0,"Channel [-]")
sheet1.write(1,1,"Lineal Energy [keV/um]")
sheet1.write(1,2,"Counts [-]")
sheet1.write(0,4,"Low Gain Spectrum")
sheet1.write(1,4,"Channel [-]")
sheet1.write(1,5,"Lineal Energy [keV/um]")
sheet1.write(1,6,"Counts [-]")
sheet1.write(0,8,"Combined Spectrum")
sheet1.write(1,8,"Lineal Energy [keV/um]")
sheet1.write(1,9,"Counts [-]")
for i in range(len(highChannel)):
    sheet1.write(i+2,0,highChannel[i])
    sheet1.write(i+2,1,highy[i])
    sheet1.write(i+2,2,highGain[i])
for i in range(len(lowChannel)):
    sheet1.write(i+2,4,lowChannel[i])
    sheet1.write(i+2,5,lowy[i])
    sheet1.write(i+2,6,lowGain[i])
for i in range(len(linBins)):
    sheet1.write(i+2,8,linBins[i])
    sheet1.write(i+2,9,linCount[i])

sheet2.write(0,0,"Lineal energy [keV/um]")
#sheet2.write(0,1,"f(y)")
sheet2.write(0,1,"yf(y)")
for i in range(len(center)):
    sheet2.write(i+1,0,center[i])
#    sheet2.write(i+1,1,fy[0][i])
    sheet2.write(i+1,1,y_fy[0][i])

sheet3.write(0,0,"Lineal energy [keV/um]")
sheet3.write(0,1,"yd(y)")
for i in range(len(center)):
    sheet3.write(i+1,0,center[i])
#   sheet3.write(i+1,1,dy[0][i])
    sheet3.write(i+1,1,y_dy[0][i])

sheet4.write(0,0,"Lineal energy [keV/um]")
sheet4.write(0,1,"y*q(y) [-]")

sheet5.write(0,0,"Cf calibration - proton edge")
sheet5.write(2,0,"Fermi parameters:")
sheet5.write(2,2,a)
sheet5.write(2,3,b)
sheet5.write(2,4,c)
sheet5.write(3,0,"inflex point:")
sheet5.write(3,1,round(c, 2))
sheet5.write(4,0,"interception point:")
sheet5.write(4,1,round(c+2/b, 2))
sheet5.write(5,0,"Proposed change of voltage:")
sheet5.write(5,2,round(CfpropVolChange, 0))
sheet5.write(7,0,"Channel [-]")
sheet5.write(7,1,"Counts*Channel^2 [-]")
for i in range(len(CfCalChn2)):
    sheet5.write(i+8,0,CfCalChn2[i])
    sheet5.write(i+8,1,CfCalCoun2[i])

wb.save(xlsName)

fig0 = plt.figure()
ax1 = fig0.add_subplot(1,1,1)
#line0 = ax1.plot(CfCalChn,CfCalCoun, 'b+', color='red', lw=2)
line0 = ax1.plot(CfCalChn,CfCalCoun, '-', color='red', lw=2)
line2 = ax1.plot(np.linspace(min(CfCalChn2)+25,max(CfCalChn2)-25,500),Fermi(np.linspace(min(CfCalChn2),max(CfCalChn2),500),a,b,c), '-', color='green', lw=2)
line3 = plt.axvline(x=c+2.5, linestyle='--', label='Finded inflex point = '+str(round(c+2.5, 1)))
#line4 = plt.axvline(x=100, linestyle=':', label='Desired position = 100.0')
ax1.set_xscale('log')
ax1.set_title(name+" - proton edge calibration")
ax1.set_ylabel("Counts*Channel^2 [-]")
ax1.set_xlabel("Channel [-]")
ax1.legend()
figName=name+"-p-edge-calibration.png"

plt.savefig(figName)

if CfpropVolChange!=0:
    print('Consider to change of voltage by approximately '+str(round(CfpropVolChange, 0))+' Volt')
if (CfpropVolChange==0):
    print('Calibration OK.')
print('Press any key to continue.')
input()

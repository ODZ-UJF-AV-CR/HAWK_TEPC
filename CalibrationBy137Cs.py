import numpy as np
import matplotlib.pyplot as plt
import math
import os
from xlwt import Workbook
from scipy.optimize import curve_fit

simulited_site_diameter = 2 # [μm]
SMScalKoef = 0.30093 # linear koeficient of SMS calibration curve
DPEcalKoef = 0.0526307 # linear koeficient of DPE calibration curve

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

# SMS evaluvation
def SMS(highChannel, highGain):
    suma1,suma2,SMS = 0, 0, 0
    for i in range(6,150):
        suma1 = suma1 + highGain[i]*highChannel[i]*highChannel[i]*highChannel[i]
        suma2 = suma2 + highGain[i]*highChannel[i]*highChannel[i]
    SMS = suma1/suma2
    return SMS

# select region of e-edge
def select_e_edge(center,y_dy,SMS):
    CsCalC, CsCalD = [],[]
    zerNm = 0
    for i in range(len(center)):
        if (center[i] >= ((SMS-1)*0.1)) and (zerNm < 10):
            CsCalC.append(center[i])
            CsCalD.append(y_dy[i])
            if (y_dy[i] == 0):
                zerNm = zerNm+1
    return list(CsCalC),list(CsCalD)

def Fermi(x, a, b, c):
    return a/(1+np.exp(b*(x-c)))

# find parameters of Fermi fit and desired position of e-edge (DPE)
def find_Fermi(CsCalC,CsCalD,simulited_site_diameter):
    def Fermi(x, a, b, c):
        return a/(1+np.exp(b*(x-c)))
    CsCalPar = curve_fit(Fermi, CsCalC, CsCalD)
    [a, b, c] = CsCalPar[0]
    DPE = 10.4*(simulited_site_diameter)**(-0.29)
    return [a, b, c], DPE

# propose approximate value of voltage change
def propose_voltage_change(SMScalKoef,SMS,DPEcalKoef,c,DPE):
    SMSpropVolChange, DPEpropVolChange = 0,0
    if abs(SMS-45.5) > 3:
        SMSpropVolChange = (SMS-45.5)/(SMScalKoef)
    if abs(c-DPE) > 0.5:
        DPEpropVolChang = (c-DPE)/(DPEcalKoef)
    return SMSpropVolChange, DPEpropVolChang

SMS = SMS(highChannel, highGain)
CsCalC, CsCalD = select_e_edge(center,y_dy[0],SMS)
[a, b, c], DPE = find_Fermi(CsCalC,CsCalD,simulited_site_diameter)
SMSpropVolChange, DPEpropVolChang = propose_voltage_change(SMScalKoef,SMS,DPEcalKoef,c,DPE)

# Ulozeni do excelovskych souboru
sheet1=wb.add_sheet("Counts")
sheet2=wb.add_sheet("Frequency distribution")
sheet3=wb.add_sheet("Dose distribution")
sheet4=wb.add_sheet("Equivalent dose distribution")
sheet5=wb.add_sheet("Calibrations")

sheet1.write(0,0,"High Gain Spectrum")
sheet1.write(1,0,"Channel [-]")
sheet1.write(1,1,"Lineal Energy [keV/μm]")
sheet1.write(1,2,"Counts [-]")
sheet1.write(0,4,"Low Gain Spectrum")
sheet1.write(1,4,"Channel [-]")
sheet1.write(1,5,"Lineal Energy [keV/μm]")
sheet1.write(1,6,"Counts [-]")
sheet1.write(0,8,"Combined Spectrum")
sheet1.write(1,8,"Lineal Energy [keV/μm]")
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

sheet2.write(0,0,"Lineal energy [keV/μm]")
sheet2.write(0,1,"yf(y)")
for i in range(len(center)): 
    sheet2.write(i+1,0,center[i])
    sheet2.write(i+1,1,y_fy[0][i])

sheet3.write(0,0,"Lineal energy [keV/μm]")
sheet3.write(0,1,"yd(y)")    
for i in range(len(center)):
    sheet3.write(i+1,0,center[i])
    sheet3.write(i+1,1,y_dy[0][i])

sheet4.write(0,0,"Lineal energy [keV/μm]")
sheet4.write(0,1,"yq(y) [-]")
for i in range(len(center)):
    sheet4.write(i+1,0,center[i])
    sheet4.write(i+1,1,y_dy[0][i])

sheet5.write(0,0,"Calibration of fitting by Fermi function")
sheet5.write(1,0,"Fermi fit parameters:")
sheet5.write(1,2,a)
sheet5.write(1,3,b)
sheet5.write(1,4,c)
sheet5.write(1,6,"SMS:")
sheet5.write(0,6,"Calibration by SMS")
sheet5.write(1,7,SMS)
sheet5.write(2,0,"Desired values:")
sheet5.write(2,4,DPE)
sheet5.write(2,7,45.5)
sheet5.write(3,0,"Proposed change of voltage:")
sheet5.write(3,7,round(SMSpropVolChange, 0))
sheet5.write(3,4,round(DPEpropVolChang, 0))
sheet5.write(5,0,"Lineal energy [-]")
sheet5.write(5,1,"yd(y) [-]")
for i in range(len(CsCalC)):
    sheet5.write(i+6,0,CsCalC[i])
    sheet5.write(i+6,1,CsCalD[i])
    
wb.save(xlsName)

fig0 = plt.figure()
ax1 = fig0.add_subplot(1,1,1)
line1 = ax1.plot(center,y_dy[0], linestyle='-', marker='+', color='red', lw=2)
line2 = ax1.plot(np.linspace(min(CsCalC),max(CsCalC),500),Fermi(np.linspace(min(CsCalC),max(CsCalC),500),a,b,c), '-', color='green', lw=2)
line3 = plt.axvline(x=c, linestyle='--', label='Finded e-edge = '+str(round(c, 2)))
line4 = plt.axvline(x=DPE, linestyle=':', label='Desired e-edge = '+str(round(DPE, 2)))
ax1.set_xscale('log')
ax1.set_title(name+" - electron edge calibration")
ax1.set_ylabel("yd(y)")
ax1.set_xlabel("y [keV/μm]")
ax1.legend()
figName=name+"-calibration.png"
plt.savefig(figName)

if SMSpropVolChange!=0:
    print('Consider to change of voltage by approximately '+str(round(SMSpropVolChange, 0))+' Volt')
if (SMSpropVolChange==0):
    print('SMS Calibration OK.')
#if DPEpropVolChang!=0:
#    print('Consider to change of voltage by approximately '+str(round(DPEpropVolChang, 0))+' Volt')
#if (DPEpropVolChang==0):
#    print('DPE Calibration OK.')
print('Press any key to continue.')
input()

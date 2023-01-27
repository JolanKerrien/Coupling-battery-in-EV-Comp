# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 09:54:18 2023

@author: Jolan

This program aim to compare a classical heating using a vapor compression cycle
with ambiant air to the same one but using the heat provided by the battery.
"""

# import pandas as pd
from pina import PinchAnalyzer, make_stream
import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP
import CoolProp.Plots as CPP
import os
script_dir = os.path.dirname(__file__)
import sys
sys.path.append(script_dir)
import Project_Ad as Ad


#######################################################################################
###############################   INITIALISATION     ##################################
#######################################################################################
effComp = 0.62
m_ref = 0.03 #kg/s
pinch = 5 #K
DeltaT_air = 5 #K #Init the DeltaTair at the HEX that not provide the cabin
T_amb = np.array([-20,-10,0,10,20,30,40])+273.15 #K
T_blown = np.array([45,39,34,27.5,20,7,5])+273.15 #K
rho_air = 1.300 #kg/m3
m_air = np.array([450,400,350,300,400,450,450])*rho_air/3600 #kg/s
m_air2 = 0.31 #kg/s #Init the m_air at the HEX that not provide the cabin
#0.31 max: 1.6
Cp_air = 1005 #J/kg/K
Q_blown = np.array([m*Cp_air*(T_blown[i]-T_amb[i]) for i,m in enumerate(m_air)]) #W
#Q_blown[:4] = Q_cond ; Q_blown[5:] = Q_evap
Cp_glyc = 3195 #J/kg/K (0.5*Cp_eau[4190]+0.5*Cp_Ethglycol[2200])

#######################################################################################
###################################   FUNCTIONS     ###################################
#######################################################################################

#FUNCTIONS PLOT P-H
def connectpoints(x,y,p1,p2): #Print the P-h diagram
    x1, x2 = x[p1], x[p2]
    y1, y2 = y[p1], y[p2]
    plt.plot([x1,x2],[y1,y2],'k-')
    
def PlotCycle(idx):
#return the P-h diagram of the cycle for a given T_amb (idx)
    plot = CPP.PropertyPlot('HEOS::R1234yf', 'PH', unit_system='EUR', tp_limits='ACHP');
    plt.ylim((1,HP[idx]*1e-5+10));
    plt.xlim((150,H2[idx]*1e-3+50));
    plot.calc_isolines(CP.iQ, num=10);
    plot.calc_isolines(CP.iT, num=10);
    PTP = np.array([BP[idx]*1e-5,HP[idx]*1e-5,HP[idx]*1e-5,BP[idx]*1e-5]) #bar
    HTP = np.array([H1[idx],H2[idx],H3_4[idx],H3_4[idx]])*1e-3 # kJ/kg
    for i in [0,1,2]:
        connectpoints(HTP, PTP, i, i+1)
    connectpoints(HTP,PTP,3,0)
    plt.scatter(HTP,PTP,s=40,c='blue');
    Mod = 'OFF'
    if idx<4:
        Mod = 'Heater'
    if idx>4:
        Mod = 'Chiller'
    plt.title(f"P-h diagram for a {Mod} using R-1234yf at {T_amb[idx]-273.15}°C")
    plot.show();



#FUNCTION COMPOSITE CURVE
def Composite(idx,pinch = pinch): #0,1,2,3,5,6
#return the composite curve of evaporator and condenser side for a given T_amb (idx)
    # Refrigeration side
    cold_Ev = make_stream(Q_evap[idx], T1_4[idx]-273.15, T1_4[idx]-273.15 ) #Evaporator
    Hsat = CP.PropsSI('H','P',HP[idx],'Q',1,'R1234yf')
    Qsens = m_ref[idx]*(H2[idx]-Hsat)
    Qlat = (Q_cond[idx])*1e3-Qsens
    hot_Condlat = make_stream((Qlat*1e-3), T3[idx]-273.15,T3[idx]-273.15)
    hot_Condsens = make_stream(Qsens*1e-3, T3[idx]-273.15, T2[idx]-273.15)
    
    #Air side
    if idx<4: 
        Mod = 'Heater'
        hot_EvAir = make_stream(-Q_evap[idx], T_amb[idx]-273.15,T_air_out[idx]-273.15)
        cold_CondAir = make_stream(-Q_cond[idx], T_amb[idx]-273.15,T_blown[idx]-273.15)
        
    elif idx>4: 
        Mod = 'Chiller'
        hot_EvAir = make_stream(-Q_evap[idx], T_amb[idx]-273.15,T_blown[idx]-273.15) 
        cold_CondAir = make_stream(-Q_cond[idx], T_amb[idx]-273.15,T_air_out[idx-1]-273.15)
    else: return 'VPC OFF'
    #Plot
    analyzerEv,analyzerCond = PinchAnalyzer(pinch/2),PinchAnalyzer(pinch/2)
    analyzerEv.add_streams(cold_Ev, hot_EvAir) #Evaporator
    analyzerCond.add_streams(cold_CondAir, hot_Condlat,hot_Condsens) #Condenser
    #Plot Evaporator
    plt.plot(*analyzerEv.hot_composite_curve, color="tab:red", linestyle="--", label="Air")
    plt.plot(*analyzerEv.cold_composite_curve, color="tab:blue", linestyle="-", label="Evaporator")
    plt.legend()
    plt.title(f"Composite curves for evaporator in {Mod} mode T_amb = {T_amb[idx]-273.15}°C")
    plt.xlabel("Heat flow [kW]")
    plt.ylabel("Actual temperature [\u2103]")
    plt.show()
    #Plot Condenser
    plt.plot(*analyzerCond.hot_composite_curve, color="tab:red", linestyle="--", label="Condenser")
    plt.plot(*analyzerCond.cold_composite_curve, color="tab:blue", linestyle="-", label="Air")
    plt.legend()
    plt.title(f"Composite curves for condenser in {Mod} mode T_amb = {T_amb[idx]-273.15}°C")
    plt.xlabel("Heat flow [kW]")
    plt.ylabel("Actual temperature [\u2103]")
    plt.show()


#FUNCTION ADJUST M_REF
def bisect_m_refHeater(j,eps,max_iter,m_ref_min,m_ref_max):
    i=0
    while i<max_iter:
        m_ref_mid = (m_ref_min+m_ref_max)/2
        q_cond = m_ref_mid*(H2[j]-H3_4[j])
        if abs(q_cond - Q_blown[j]) <eps:
            return m_ref_mid
        elif q_cond-Q_blown[j]>eps:
            m_ref_max = m_ref_mid
        else:
            m_ref_min = m_ref_mid
        i+=1
    return m_ref_mid

def bisect_m_refChiller(j,eps,max_iter,m_ref_min,m_ref_max):
    i=0
    while i<max_iter:
        m_ref_mid = (m_ref_min+m_ref_max)/2
        q_ev = m_ref_mid*(H3_4[j]-H1[j])
        if abs(q_ev - Q_blown[j]) <eps:
            return m_ref_mid
        elif q_ev-Q_blown[j]>eps:
            m_ref_min = m_ref_mid
        else:
            m_ref_max = m_ref_mid
        i+=1
    return m_ref_mid

def Adjustm_ref(eps = 1e-2,max_iter=1000,m_ref_min = 0.01,m_ref_max = 0.1):
    m_ref = []
    for idx in range(7):
        if idx<=3:
            m_ref.append(bisect_m_refHeater(idx, eps,max_iter,m_ref_min,m_ref_max))
        elif idx ==4:
            m_ref.append(0)
        else:
            m_ref.append(bisect_m_refChiller(idx, eps,max_iter,m_ref_min,m_ref_max))
    return m_ref


#FUNCTIONS ADJUST M_AIR2 
def bisect_m_air2Heat(j, eps,max_iter,m_air_min,m_air_max,pinch):
    i = 0
    while i < max_iter:
        i += 1
        m_air_mid = (m_air_min + m_air_max) / 2
        T_air_out = T_amb[j] + Q_evap[j]*1e3 / m_air_mid / Cp_air
        if abs(T_air_out - T1_4[j]-pinch) < eps and (T_air_out - T1_4[j]-pinch)>0: 
            return m_air_mid
        elif T_air_out - T1_4[j]-pinch > eps:
            m_air_max = m_air_mid
        else:
            m_air_min = m_air_mid
    return m_air_mid

def bisect_m_air2Cold(j, eps,max_iter,m_air_min,m_air_max,pinch):
    i = 0
    while i < max_iter:
        i += 1
        m_air_mid = (m_air_min + m_air_max) / 2
        T_air_out = T_amb[j+1] + Q_sens[j+1]*1e3 / m_air_mid / Cp_air #pinch at saturation
        if abs(T3[j+1]-T_air_out - pinch) < eps and T3[j+1]-T_air_out - pinch>0: #so that it a bit more than pinch
            return m_air_mid
        elif T3[j+1]-T_air_out -pinch > eps:
            m_air_max = m_air_mid
        else:
            m_air_min = m_air_mid
    return m_air_mid

def Adjustm_air2(eps = 1e-2,max_iter=1000,m_air_min = 0.01,m_air_max = 2, pinch = pinch):
    m_air2 = []
    for idx in range(6):
        if idx<=3: #Heating Mod --> Looking at the cold side(evap)
            m_air2.append(bisect_m_air2Heat(idx, eps,max_iter,m_air_min,m_air_max,pinch))
        else: # Chiller mod --> Looking at the hot side (cond)
            m_air2.append(bisect_m_air2Cold(idx, eps,max_iter,m_air_min,m_air_max,pinch))
    return m_air2

#######################################################################################
###############################   Basic requirements    ###############################
#######################################################################################
# Point 3
T3 = np.concatenate((T_blown[:4],T_amb[4:]+DeltaT_air))+pinch #K
HP = CP.PropsSI('P','T',T3,'Q',0,'R1234yf') #Pa
H3_4 = CP.PropsSI('H','T',T3,'Q',0,'R1234yf') #J/kg
S3 = CP.PropsSI('S','T',T3,'Q',0,'R1234yf') #J/kg/K

# Point 1
T1_4 = np.concatenate((T_amb[:4]-DeltaT_air,T_blown[4:]))-pinch #K
BP = CP.PropsSI('P','T',T1_4,'Q',1,'R1234yf') #Pa
S1 = CP.PropsSI('S','T',T1_4,'Q',1,'R1234yf') #J/kg/K
H1 = CP.PropsSI('H','T',T1_4,'Q',1,'R1234yf') #J/kg

# Point 4
#T4 = T1
#H4 = H3
x4 = CP.PropsSI('Q',"H",H3_4,'P',BP,'R1234yf')
S4 = CP.PropsSI('S','T',T1_4,'Q',x4,'R1234yf') #J/kg/K

# Point 2
H2isen = CP.PropsSI("H","S",S1,"P",HP,'R1234yf') #J/kg
H2 = (H2isen -H1*(1-effComp))/effComp #J/kg
T2 = CP.PropsSI('T','H',H2,'P',HP,'R1234yf') #K
S2 = CP.PropsSI('S','T',T2,'P',HP,'R1234yf') #J/kg/K

#Plot the entire cycle
# PlotCycle(5) #Plot whatever you want: 0,1,2,3 is heating mod and 5,6 chiller mod

#######################################################################################
################################   Energy Balance     #################################
#######################################################################################
m_ref = Adjustm_ref()
# print(f'm_ref = {m_ref} kg/s')
Q_cond = m_ref*(H2-H3_4)*1e-3 #kW
Q_evap = m_ref*(H3_4-H1)*1e-3 #kW
W = m_ref*(H2-H1)*1e-3 #kW
COP = np.concatenate((Q_cond[:4]/W[:4],((-1)*Q_evap[5:]/W[5:])))

#Calcul Hsat to adjust m_air
H_sat = CP.PropsSI('H','P',HP,'Q',1,'R1234yf') #J/kg
Q_sens = m_ref*(H_sat-H3_4)*1e-3 #kW

#Adjust m_air2
m_air2 = Adjustm_air2()
T_air_out = np.concatenate((T_amb[:4]+Q_evap[:4]*1e3/m_air2[:4]/Cp_air,
                                T_amb[5:]+Q_cond[5:]*1e3/m_air2[4:]/Cp_air))
# print(f'm_air2 = {m_air2} kg/s')

#Plot Graph vs Tamb
# plt.plot(T_amb-273.15,np.insert(COP,4,0),'.') #Vapor Compression cycle off
# plt.axvline(x=20,color='gray',linestyle='--') 
# plt.xlabel("T_amb [°C]")
# plt.ylabel("COP")
# plt.savefig(os.path.join(script_dir, 'COP vs Tamb.png'),dpi=200);
# plt.show()
# plt.plot(T_amb-273.15,W,'.',color = 'orange')
# # plt.plot(T_amb-273.15,abs(Q_blown)*1e-3,'.')
# plt.axvline(x=20,color='gray',linestyle='--') 
# plt.xlabel("T_amb [°C]")
# plt.ylabel("W_comp [kW]")
# plt.savefig(os.path.join(script_dir, 'W_comp vs Tamb.png'),dpi=200);
# plt.show()

#######################################################################################
###############################   Composite curve     #################################
#######################################################################################

# Composite(5) #Plot whatever you want: 0,1,2,3 is heating mod and 5,6 chiller mod

#######################################################################################
###############################   Advanced requirements     ###########################
#######################################################################################

#import Advanced Requirements outputs
COPAd,WAd,Q_condAd, Q_evapAd, m_refAd,m_glycAd = Ad.Advanced()

#Plot COP vs T_amb
plt.plot(T_amb-273.15,np.insert(COP,4,0),'.',color = 'red',label = 'HVAC')
plt.plot(T_amb-273.15,np.insert(COPAd,4,0),'.', color ='blue',label = 'Battery coupled HVAC')
plt.axvline(x=20,color='gray',linestyle='--')
plt.title('COP vs ambient temperature for 2 scenarios')
plt.xlabel("T_amb [°C]")
plt.ylabel("COP")
plt.legend()
# plt.savefig(os.path.join(script_dir, 'COP vs Tamb.png'),dpi=200);
plt.show()

#Plot W vs T_amb
plt.plot(T_amb-273.15,W,'.',color = 'red', label = 'HVAC')
plt.plot(T_amb-273.15,WAd,'.',color = 'blue',label = 'Battery coupled HVAC')
# plt.plot(T_amb-273.15,abs(Q_blown)*1e-3,'.')
plt.axvline(x=20,color='gray',linestyle='--')
plt.title('Work vs ambient temperature for 2 scenarios') 
plt.xlabel("T_amb [°C]")
plt.ylabel("W_comp [kW]")
plt.legend()
# plt.savefig(os.path.join(script_dir, 'W_comp vs Tamb.png'),dpi=200);
plt.show()

#Plot Q vs T_amb
plt.plot(T_amb-273.15,Q_cond,'.',color = 'darkred', label = 'HVAC (cond)')
plt.plot(T_amb-273.15,Q_condAd,'.',color = 'darkblue',label = 'Battery coupled HVAC (cond)')
plt.plot(T_amb-273.15,Q_evap,'.',color = 'tomato', label = 'HVAC (evap)')
plt.plot(T_amb-273.15,Q_evapAd,'.',color = 'cyan',label = 'Battery coupled HVAC (evap)')
# plt.plot(T_amb-273.15,abs(Q_blown)*1e-3,'.')
plt.axvline(x=20,color='gray',linestyle='--')
plt.title('Heat vs ambient temperature for 2 scenarios') 
plt.xlabel("T_amb [°C]")
plt.ylabel("Q [kW]")
plt.legend()
plt.show()

#Plot ref flow rate vs T_amb:
plt.plot(T_amb-273.15,m_ref,'.',color = 'red',label = 'HVAC')
plt.plot(T_amb-273.15,m_refAd,'.',color = 'blue', label = 'Battery coupled HVAC')
plt.axvline(x=20,color='gray',linestyle='--') 
plt.title('refrigerant mass flow rate vs ambient temperature for 2 scenarios')
plt.xlabel("T_amb [°C]")
plt.ylabel("mass flow rate [kg/s]")
plt.legend()
# plt.savefig(os.path.join(script_dir, 'm_ref vs Tamb.png'),dpi=200);
plt.show()

#Plot Heat capacity rate vs T_amb
plt.plot(T_amb-273.15,np.insert(Cp_air*np.array(m_air2),4,0),'.',color = 'red',label = 'HVAC (air)')
plt.plot(T_amb-273.15,np.insert(Cp_glyc*m_glycAd,4,0),'.',color = 'blue', label = 'Battery coupled HVAC (glyc)')
plt.axvline(x=20,color='gray',linestyle='--')
plt.title('Heat capacity rate vs ambient temperature for 2 scenarios')
plt.xlabel("T_amb [°C]")
plt.ylabel("Cp*m [W/K]")
plt.legend()
# plt.savefig(os.path.join(script_dir, 'mCp vs Tamb.png'),dpi=200);
plt.show()













    
    
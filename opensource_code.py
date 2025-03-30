from numpy import pi
import numpy as np
from scipy.integrate import quad
from scipy.optimize import fsolve
import math
import sys
from scipy.optimize import root
import matplotlib.pyplot as plt
# from radconv1d import make_profiles
# from radconv1d import radconv
from openpyxl import Workbook
import pandas as pd

# >> input parameters
T_s = float(input("T_s (in K) = ")) # Temperature of the surface of magma ocean in K (e.g., 7154)
T_mb = float(input("T_mb (in K) = ")) # Temperature of the bottom of magma ocean's boundary layer in K (e.g., 7354)
T_ab = float(input("T_ab (in K) = ")) # Temperature of the top of atmosphere's boundary layer in K (e.g., 7154)
R_mo = float(input("R_mo (in km) = ")) # Depth of the magma ocean in km (e.g., 2000)
M = float(input("M (in kg) = ")) # Mass of the proto-Earth in e24 kg (e.g., 5.9202*1e24 )
R0 = float(input("R0 (in km) = ")) # Radius of the proto-Earth in km (e.g., 6371)
Ma = float(input("Ma (in kg) = ")) # Mass of the primary atmosphere in e21kg (e.g., 5.9202* 1e21)
# << input parameters

k_mb = 2  # Thermal conductivity in W*K-1*m-1
rou_mb = 4000  # Magma ocean density in kg*m-3
alpha_mb = 5e-5  # Volume expansion coefficient in K-1
eta_mb = 0.1  # Dynamic viscosity in pa*s 
cp_mb = 5000  # Specific heat of the mantle in J*K*kg-1
G = 6.67259e-11  # Gravitational constant in m3*kg-1*s-2
K = 1  # Reduction factor in m2*s-1
kB = 1.381e-23  # Boltzmann's constant in m2*kg*s-2*K-1
gamma = 5 / 3  # Heat capacity ratio
gamma_0 = 4 / 3  # Reference heat capacity ratio
NA = 6.02214076 * 10 ** 23  # Avogadro's constant
R = 8.314  # Gas constant in J* mol-1 *K-1
xi1 = 0.5  # The efficiency of mass loss
kappa_th = 1  # Thermal opacity of molecular hydrogen in m2*kg-1
sigma = 5.67 * 10 ** (-8)  # Stefan-Boltamann constant in W*m-2*K-4
aE = 1e3 * 1.496 * 1e8  # Semi-major axis of the Earth in m
M_sun = 1.989e30  # Mass of the Sun in kg
L = 0.7 * 3.845 * 1e26  # Luminosity of the young Sun

delta_t = 1  # Time step in Myrs

# >>> molecular weight in g
m_H = 1.008 * 1.993 / 12 * 10 ** (-23)  
m_O = 16 * 1.993 / 12 * 10 ** (-23)
m_He = 4.0026 * 1.993 / 12 * 10 ** (-23)
Mol_CO2 = 44.010
m_H2O = 2 * m_H + m_O
m_CO2 = Mol_CO2 * 1.993 / 12 * 10 ** (-23)
m_20Ne = 19.99 * 1.993 / 12 * 10 ** (-23)
m_22Ne = 21.99 * 1.993 / 12 * 10 ** (-23)
m_36Ar = 35.97 * 1.993 / 12 * 10 ** (-23)
m_38Ar = 37.96 * 1.993 / 12 * 10 ** (-23)
# <<< molecular weight in g

g = (G * M) / ((R0 * 1000) ** 2)
M_mo = rou_mb * 4 * pi * ((R0 * 1000) ** 3 - ((R0 - R_mo) * 1000) ** 3) / 3 
initial_H = Ma / M
Teq = (L / (16 * pi * sigma * (aE ** 2))) ** (1 / 4) 
MO_Amount_Ar_36 = 0
MO_Amount_Ar_38 = 0
MO_Amount_He = 0
MO_Amount_Ne_20 = 0
MO_Amount_Ne_22 = 0
M_H2O = 0 
M_CO2 = 0 
M_H2 = 1 * Ma
H_1bar = (1 * 1e5) / g * 4 * pi * (R0 * 1e3) ** 2 * 1e3 / m_H
Amount_O = M_H2O * 1e3 / (2 * m_H + m_O)
Amount_H_nebula = initial_H * M * 1e3 / m_H  
Amount_H = 2 * Amount_O + Amount_H_nebula
Amount_CO2 = M_CO2 * 1e3 / m_CO2
Amount_Ne_20 = 6.46e-4 * Amount_H_nebula + MO_Amount_Ne_20
Amount_Ar_36 = 1.31e-5 * Amount_H_nebula + MO_Amount_Ar_36
Amount_He = 9.88e-2 * Amount_H_nebula + MO_Amount_He
Amount_Ne_22 = Amount_Ne_20 / 13.75 + MO_Amount_Ne_22
Amount_Ar_38 = Amount_Ar_36 / 5.50005 + MO_Amount_Ar_38
g_mb = g 
g_ab = g 

Fs = 0.089 * ((k_mb * rou_mb * (T_s - T_mb) ** 2) ** (2 / 3)) * (cp_mb * g_mb * alpha_mb / eta_mb) ** (
            1 / 3)  # Thermal flux of magma ocean in kg*s-3
R_acr = 1000  # Critical Rayleigh number
alpha_ab = 3.66e-3  # Volume expansion parameter of H2 
cp_ab = 1.43e4  # Specific heat capacity at constant pressure of H2 in J*kg-1*K-1

d = 2*5.3e-11 # Mean free path
miu = ((Amount_O * m_O + Amount_H * m_H + Amount_He * m_He
        + Amount_CO2 * m_CO2 + Amount_Ne_20 * m_20Ne + Amount_Ne_22 * m_22Ne + Amount_Ar_36 * m_36Ar + Amount_Ar_38 * m_38Ar)
    / (Amount_O + (
            Amount_H) + Amount_CO2 + Amount_Ne_20 + Amount_Ne_22 + Amount_Ar_36 + Amount_Ar_38 + Amount_He))  
miu0 = miu
rou_ab = (7.2425E27/T_s)*miu*1e-3 * (Ma*g/ (4 * pi * (R0* 1e3) ** 2) * 1e-5)
gamma = 5 / 3

y1 = 1.8 * (3 / 2)  
y2 = (d / (1e-10)) ** (-2)
y3 = (T_ab / 1000) ** (1 / 2)
y4 = (miu * 1e-3 / (1.661 * (1e-27))) ** (-1 / 2)
k_ab = y1 * y2 * y3 * y4
eta_ab = 8.4 * (1e-5) * (d / (1e-10)) ** (-2) * (T_ab / 1000) ** (1 / 2) * (miu * (1e-3) / (1.661 * (1e-27))) ** (
            1 / 2)  
T_ab = T_s - ((R_acr * eta_ab * (Fs) ** 3) / (cp_ab * (rou_ab ** 2) * g_ab * alpha_ab * (k_ab ** 2))) ** (
            1 / 4)  
R_ab = R0 * 1000 + ((R_acr * eta_ab * (k_ab ** 2)) / (cp_ab * (Fs) * g_ab * alpha_ab * (rou_ab ** 2))) ** (1 / 4)  # m
T_B = 2 / (3 - gamma) * (T_ab - ((gamma - 1) * G * M * (miu0 * (1e-3))) / ((gamma * R_ab * kB)))  # Bondi radius' temperature
T_B1 = 48 * ((gamma / gamma_0) ** (-1)) * (miu0 * (1e-3) / (1.661 * (1e-27)))  
T_B2 = ((2 + 3 ** (1 / 2)) / 2) ** (1 / 4) * Teq
R_B = 2 * G * M * (miu0 * (1e-3)) / (gamma * kB * T_B)  
T_abc1 = (Teq * (3 - gamma) * ((2 + (3) ** (1 / 2)) / 2) ** (1 / 4)) / 2 + ((gamma - 1) * G * M * (miu0 * (1e-3))) / (
            gamma * kB * R_ab)  # Bondi radius' critical temperature in K, if T_B<T_BC
T_abc2 = 2 * G * M * miu0 * 1e-3 / (gamma * kB * R_ab)
print("T_abc1 = ", T_abc1)
print("T_abc2 = ", T_abc2)
rou_B = rou_ab * np.exp(-G * M * (miu0 * 1e-3) * (1 / R_ab - 1 / R_B) / (kB * T_B))
tau_B = kappa_th * rou_B * kB * T_B * (R_B) ** 2 / (G * M * (miu0 * 1e-3))
# tau_B = 2*kappa_th*rou_B*R_B/gamma
FB = 2 * sigma * gamma * (T_B ** 4 - (2 + 3 ** (1 / 2)) / 4 * Teq ** 4) / (3 * kappa_th * rou_B * R_B + gamma)
T_rcb = ((2 + 3 ** (1 / 2)) / 2) ** (1 / 4) * Teq
Frcb = 4 * sigma * (T_rcb ** 4 - (2 + 3 ** (1 / 2)) / 4 * Teq ** 4) / (3 * tau_B + 2)

Ma_c = 8 * pi * (R0 * 1e3) ** 4 / (3 * kappa_th * R_B ** 2) * np.exp(
    gamma * (R_B / (R0 * 1e3) - 1) / 2)  # The lowest atmospheric mass of regime 1

def integrand(r):
    return ((((3 - gamma) / 2) + ((gamma - 1) * R_B) / (2 * r)) ** (1 / (gamma - 1)) * (r ** 2))

Ma_c1 = 4 * pi * gamma / (3 * kappa_th * R_B) * quad(integrand, R_ab, R_B)[0]

Rx = 1*R_B
u_e = (2*G*M/Rx)**(1/2)
print("u_e = ", u_e)
Tx = miu0 * 1e-3 * u_e ** 2 / (3 * kB)
print("Tx = ", Tx)
if Tx <= Teq:
    Tx = Teq
rou_e = G*M*(miu0*1E-3)**2/(2**(1/2) * pi * d**2 * kB * Tx * Rx**2)
Rx = (G*M*(miu0*1e-3)**2 / (2**(1/2) * pi*d**2*rou_e*kB*Tx))**(1/2)
print("Tx = ", Tx)
ratio_R = R_B / Rx
print("*ratio_R = ", ratio_R)

h = 0
j = 0
iR = 0
MA = [0] * 1000000  
j_fig = []
Ma_fig = []
T_ab_fig = []
T_B_fig = []
T_mo_fig = []
ratio_R_fig = []
Data = []

while j < 3:
    MA[j] = Ma
    j = j + 1
j = 0

if T_ab > T_abc1 and Ma > Ma_c and T_ab < T_abc2:

    while MA[j + 1] > max(Ma_c1, Ma_c) and R_B > R_ab and T_ab > T_abc1 and T_B > max(T_B1, T_B2):
        # regime 1

        def integrand(r):
            return ((((3 - gamma) / 2) + ((gamma - 1) * R_B) / (2 * r)) ** (1 / (gamma - 1)) * (r ** 2))

        Ma1 = rou_B * 4 * pi * quad(integrand, R_ab, R_B)[0]  
        MA[j + 1] = min(Ma1, MA[j + 1])
        if j >= 1 and MA[j + 1] != MA[j] - M_B:
            MA[j + 1] = MA[j] - M_B
        print("*Ma = ", MA[j+1])

        u_B = (gamma * kB * T_B / (miu0 * (1e-3))) ** (1 / 2)  # Sonic speed at Bondi Radius in m*s-1
        M_B = 4 * pi * xi1 * rou_B * u_B * (R_B ** 2)  # mass loss rate at Bondi radius in kg/s
        print("*MB = ", M_B)

        MA[j + 2] = MA[j + 1] - M_B
        PPa =  MA[j+1]*g*1e-5/(4*pi*(R0*1e3)**2)
        print("*Ma = ", MA[j+1])
        print("*P = ", MA[j+1]*g*1e-5/(4*pi*(R0*1e3)**2))

        F_B = 2 * sigma * gamma * (T_B ** 4 - (2 + 3 ** (1 / 2)) / 4 * Teq ** 4) / (3 * kappa_th * rou_B * R_B + gamma)
        if F_B < 0:
            F_B = 0
            print("F_B = 0")
            break
        
        Fml = G * M * M_B / (4 * pi * R_B ** 3)  # mass-loss rate

        Fs = F_B
        print("*Fs = ",Fs)
        print("*Fml/Fs = ",f"{Fml/Fs:.3e}")
        deltaT_mb = 4 * pi * ((R0 * 1e3) ** 2) * Fs / (M_mo * cp_mb)  
        print(deltaT_mb)
        T_mb = T_mb - deltaT_mb  

        T_s = T_mb - (((1 / (k_mb * rou_mb)) * (Fs / 0.089 * (cp_mb * g_mb * alpha_mb / eta_mb) ** (-1 / 3)) ** (
                    3 / 2)) ** (1 / 2)) - 200

        rou_B = MA[j + 2] / (4 * pi * quad(integrand, R_ab, R_B)[0])
        rou_ab = rou_B * ((3 - gamma) / 2 + (gamma - 1) * R_B / (2 * R_ab)) ** (1 / (gamma - 1))

        T_ab = T_s - ((R_acr * eta_ab * (Fs) ** 3) / (cp_ab * (rou_ab ** 2) * g_ab * alpha_ab * (k_ab ** 2))) ** (
                    1 / 4)  
        if T_ab < 0:
            T_ab = 0
            print("T_ab = 0")
            break
        T_ab_fig.append(T_ab)

        if j > 1 and T_ab_fig[j] > T_ab_fig[j - 1]:
            T_ab_fig[j] = min(T_ab_fig[j - 1], T_ab_fig[j])
        T_ab = T_ab_fig[j]

        R_ab = R0 * 1000 + ((R_acr * eta_ab * (k_ab ** 2)) / (cp_ab * Fs * g_ab * alpha_ab * (rou_ab ** 2))) ** (
                    1 / 4)  # in m

        k_ab = 1.8 / (gamma - 1) * ((d / (1e-10)) ** (-2)) * ((T_ab / 1000) ** (1 / 2)) * (
                    miu * (1e-3) / (1.661 * (1e-27))) ** (-1 / 2)  
        eta_ab = 8.4 * (1e-5) * ((d / (1e-10)) ** (-2)) * ((T_ab / 1000) ** (1 / 2)) * (
                    miu * (1e-3) / (1.661 * (1e-27))) ** (1 / 2)  

        T_B = 2 / (3 - gamma) * (T_ab - ((gamma - 1) * G * M * (miu0 * (1e-3))) / (
        (gamma * R_ab * kB)))  
        if T_B < 0:
            T_B = 0
            print("T_B = 0")
            break
        T_B_fig.append(T_B)
        if j > 1 and T_B_fig[j] > T_B_fig[j - 1]:
            T_B_fig[j] = min(T_B_fig[j - 1], T_B_fig[j])
        T_B = T_B_fig[j]

        R_B = 2 * G * M * (miu0 * (1e-3)) / (gamma * kB * T_B)  

        tau_B = kappa_th * rou_B * kB * T_B * (R_B) ** 2 / (G * M * (miu0 * 1e-3))

        Ma_c = 8 * pi * (R0 * 1e3) ** 4 / (3 * kappa_th * R_B ** 2) * np.exp(
            gamma * (R_B / (R0 * 1e3) - 1) / 2)  

        def integrand(r):
            return ((((3 - gamma) / 2) + ((gamma - 1) * R_B) / (2 * r)) ** (1 / (gamma - 1)) * (r ** 2))

        Ma_c1 = 4 * pi * gamma / (3 * kappa_th * R_B) * quad(integrand, R_ab, R_B)[0]

        u_e = (2 * G * M / Rx) ** (1 / 2)
        print("u_e = ", u_e)
        Tx = miu0 * 1e-3 * u_e ** 2 / (3 * kB)
        rou_e = G * M * (miu0 * 1E-3) ** 2 / (2 ** (1 / 2) * pi * d ** 2 * kB * Tx * Rx ** 2)
        Rx = (G * M * (miu0 * 1e-3) ** 2 / (2 ** (1 / 2) * pi * d ** 2 * rou_e * kB * Tx))**(1/2)
        print("Tx = ", Tx)
        print("*Rx = ", Rx)
        ratio_R = R_B / Rx
        ratio_R_fig.append(ratio_R)
        print("*ratio_R = ", ratio_R)
        T_mo_fig.append(T_s)

        if R_B > Rx:
            h = h + 1
            break

        j = j + 1
        print(j)
        j_fig.append(j)
        print("--------------------------------------")

j = j + 1
print("Time = ", j / 3600, "H")
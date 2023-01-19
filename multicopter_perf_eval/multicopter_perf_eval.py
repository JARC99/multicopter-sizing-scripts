"""Multicopter performance evaluation script"""


import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from fluids import ATMOSPHERE_1976
from tabulate import tabulate

from scipy.optimize import fsolve

sns.set_theme(style="whitegrid")
GRAPHICS_PATH = "output/"
DPI = 300

#### %% Mission Profile

endrce = 20  # min
t_charge = 3  # h

# %% Operational conditions

h = 10  # AMSL
delta_T = 10.065  # K

atm = ATMOSPHERE_1976(h, delta_T)
rho = atm.rho  # kg/m^3
p = atm.P  # Pa
T = atm.T  # K
g = atm.g  # m/s^2

# %% Aircraft Constants

n_r = 4

# Propeller constants
A_p = 5
eps_p = 0.85
lambda_p = 0.75
xi_p = 0.5
e_p = 0.83
Cfd_p = 0.015
alpha0_p = 0
K0_p = 6.11

D_p = 9.5 * 0.0254  # m - Diameter
H_p = 5 * 0.0254  # m - Pitch
B_p = 2  # Blade number

# Motor
KV0 = 920  # rpm/V - Motor KV
Im_max = 13  # A - Maximum current
Im0 = 0.5  # A - No-load current
Um0 = 10  # V - No-load voltage
Rm = 0.142  # Ohm - Motor resistance

# ESC
Ie_max = 15  # A - Maximum current
Re = 0.008  # Ohm - ESC resistance

Ictrl = 1  # Controller current

# Battery
n_cells = 4
Cbat = 2000  # np.arange(100, 10000, 10)  # 2200  # mAh - Capacity
Cmin = 0.2*Cbat  # Minimum capacity of a battery set
Rbat = 0.01  # Battery resistance
Ubat = 3.7*n_cells  # V - Total voltage
Kbat = 45*Cbat  # Maximum discharge rate
e_bat = 150  # Wh/kg

# Drag approximation
S = np.pi*D_p**2/4 * 4

C1 = 3
C2 = 1.5

# Weights
Wf_e = 0.4
W_pyld = 0.200  # kg (Camera = 0.056 kg; Gimbal = ?)
W_bat = (Cbat/1000)*Ubat/e_bat  # kg
W0 = (W_pyld + W_bat)*g/(1 - Wf_e)  # N4
# kg

MTOW = W0/g
W_e = Wf_e*MTOW


# %% Hover flight

CT = 0.25*np.pi**3*lambda_p*xi_p**2*B_p*K0_p * \
    (eps_p*np.arctan(H_p/(np.pi*D_p)) - alpha0_p)/(np.pi*A_p + K0_p)
Cd = Cfd_p + np.pi*A_p*K0_p**2/e_p * \
    (eps_p*np.arctan(H_p/(np.pi*D_p)) - alpha0_p)**2/(np.pi*A_p + K0_p)**2
CM = 1/(8*A_p)*np.pi**2*Cd*xi_p**2*lambda_p*B_p**2

T_hvr = W0/n_r
N_hvr = 60*np.sqrt(T_hvr/(rho*D_p**4*CT))
M_hvr = rho*D_p**5*CM*(N_hvr/60)**2

Um_hvr = Rm*(M_hvr*KV0*Um0/(9.55*(Um0 - Im0*Rm)) + Im0) + \
    (Um0 - Im0*Rm)/(KV0*Um0)*N_hvr
Im_hvr = M_hvr*KV0*Um0/(9.55*(Um0 - Im0*Rm)) + Im0

sigma_hvr = (Um_hvr + Im_hvr*Re) / Ubat  # [0, 1]
Ie_hvr = sigma_hvr*Im_hvr  # Missing flag
#Ie_flag = Ie_hvr > Ie_max
#     print("ESC current limit has been exceeded!")

Ibat_hvr = n_r*Ie_hvr + Ictrl
# if Ibat_hvr > Cbat*Kbat:
#     print("Battery current limit has been exceeded!")
Ue_hvr = Ubat - Ibat_hvr*Rbat

t_bat_hvr = (Cbat - Cmin)/Ibat_hvr * 60/1000  # Missing flag


# fig, (ax, ax1) = plt.subplots(2, 1, sharex=True, dpi=DPI)
# ax.plot(Cbat, t_bat_hvr)
# ax.set_xlim(left=0)
# ax.set_ylim(bottom=0)
# ax.set_ylabel("Hover Time, min")

# ax1.plot(Cbat, MTOW)
# ax1.plot(Cbat, W_e)
# ax1.plot(Cbat, W_bat)
# #ax1.axhline(max_thrust/g/T2W_rat, linestyle = "--", color="red", alpha=0.5)
# ax1.set_xlim(left=0)
# ax1.set_ylim(bottom=0)
# ax1.set_xlabel("Battery Capacity, mAh")
# ax1.set_ylabel("Mass, kg")
# ax1.legend([r"$\mathdefault{W_{0}}$",
#            r"$\mathdefault{W_{e}}$", r"$\mathdefault{W_{bat}}$"])

table_hvr = [["Endurance", round(t_bat_hvr, 3), "min"],
             ["Throttle", round(sigma_hvr, 3) * 100,
              "%"], ["ESC Current", round(Ie_hvr, 3), "A"],
             ["ESC Voltage", round(Ue_hvr, 3), "V"], [
    "Battery Current", round(Ibat_hvr, 3), "A"],
    ["Motor Angular Speed", round(N_hvr, 3), "rpm"]]

print("\nHover Flight:\n")
print(tabulate(table_hvr, headers=[
      "Parameter", "Value", "Unit"], tablefmt="github"))

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, dpi=DPI)
ax1.plot(Cbat, t_bat_hvr)
ax1.set_xlim(left=0)
ax1.set_ylim(bottom=0)
ax1.set_ylabel("Hover Time, min")

ax2.plot(Cbat, MTOW)
ax2.plot(Cbat, W_e)
ax2.plot(Cbat, W_bat)
#ax2.axhline(max_thrust/g/T2W_rat, linestyle = "--", color="red", alpha=0.5)
ax2.set_xlim(left=0)
ax2.set_ylim(bottom=0)
ax2.set_xlabel("Battery Capacity, mAh")
ax2.set_ylabel("Mass, kg")
ax2.legend([r"$\mathdefault{W_{0}}$",
            r"$\mathdefault{W_{e}}$", r"$\mathdefault{W_{bat}}$"])

# table_hvr = [["Endurance", round(t_bat_hvr, 3), "min"],
#               ["Throttle", round(sigma_hvr, 3) * 100,
#               "%"], ["ESC Current", round(Ie_hvr, 3), "A"],
#               ["ESC Voltage", round(Ue_hvr, 3), "V"], [
#     "Battery Current", round(Ibat_hvr, 3), "A"],
#     ["Motor Angular Speed", round(N_hvr, 3), "rpm"]]

# print("\nHover Flight:\n")
# print(tabulate(table_hvr, headers=[
#       "Parameter", "Value", "Unit"], tablefmt="github"))


# %% Maximum Thrust

sigma_maxT = 1.0


def eqn_system(p):
    N_maxT, M_maxT, Um_maxT, Im_maxT = p
    return (Um_maxT + Im_maxT*Re - sigma_maxT*Ubat, M_maxT - rho*D_p**5*CM*(N_maxT/60)**2,
            Um_maxT - Rm*Im0 - Rm*KV0*Um0 /
            (9.55*(Um0 - Im0*Rm))*M_maxT - ((Um0 - Im0*Rm)/(KV0*Um0))*N_maxT,
            Im_maxT - KV0*Um0/(9.55*(Um0 - Im0*Rm))*M_maxT - Im0)


N_maxT, M_maxT, Um_maxT, Im_maxT = fsolve(eqn_system, [9000, 0.12, 12, 12])

Ie_maxT = sigma_maxT*Im_maxT

Ibat_maxT = Ie_maxT*n_r + Ictrl
Ue_maxT = Ubat - Ibat_maxT*Rbat

eta_maxT = 2*np.pi/60*n_r*M_maxT*N_maxT/(Ubat*Ibat_maxT)

table_maxT = [["ESC Current", round(Ie_maxT, 3), "A"],
              ["ESC Voltage", round(Ue_maxT, 3), "V"],
              ["Battery Current", round(Ibat_maxT, 3), "A"],
              ["Motor Angular Speed", round(N_maxT, 3), "rpm"],
              ["Efficiency", round(eta_maxT, 3)*100, "%"]]

print("\nMax Thrust:\n")
print(tabulate(table_maxT, headers=[
      "Parameter", "Value", "Unit"], tablefmt="github"))

# %% Forward Flight

sigma_fwd = 0.80


def eqn_system(p):
    N_fwd, M_fwd, Um_fwd, Im_fwd = p
    return (Um_fwd + Im_fwd*Re - sigma_fwd*Ubat,
            M_fwd - rho*D_p**5*CM*(N_fwd/60)**2, Um_fwd - Rm*Im0 - Rm*KV0*Um0 /
            (9.55*(Um0 - Im0*Rm))*M_fwd - ((Um0 - Im0*Rm)/(KV0*Um0))*N_fwd,
            Im_fwd - KV0*Um0/(9.55*(Um0 - Im0*Rm))*M_fwd - Im0)


N_fwd, M_fwd, Um_fwd, Im_fwd = fsolve(eqn_system, [9000, 0.12, 12, 12])

Ie_fwd = sigma_fwd*Im_fwd
Ibat_fwd = Ie_fwd*n_r + Ictrl
Ue_fwd = Ubat - Ibat_fwd*Rbat

eta_fwd = 2*np.pi/60*n_r*M_fwd*N_fwd/(Ubat*Ibat_fwd)

T_fwd = rho*D_p**4*CT*(N_fwd/60)**2
Wpyl_max = T_fwd*n_r - W0
theta_max = np.arccos(W0/(T_fwd*n_r))

table_fwd = [["Maximum Thrust", round(T_fwd*n_r, 3), "N"],
             ["Maximum Payload", round(Wpyl_max/g, 3), "kg"],
             ["Maximum Pitch Angle", round(theta_max*180/np.pi, 3), "°"]]

print("\nMax Operational Thrust (80%):\n")
print(tabulate(table_fwd, headers=[
      "Parameter", "Value", "Unit"], tablefmt="github"))

# %% Maximum Forward Speed and Distance

theta = 30 * np.pi/180

V_maxV = np.sqrt(2*W0*np.tan(theta)/(rho*S *
                                     (C1*(1 - np.cos(theta)**3) +
                                      C2*(1 - np.sin(theta)**3))))

N_maxV = 60 * np.sqrt(W0/(rho*CT*D_p**4*n_r*np.cos(theta)))

M_maxV = W0*CM*D_p/(CT*n_r*np.cos(theta))

Um_maxV = Rm*(M_maxV*KV0*Um0/(9.55*(Um0 - Im0*Rm)) + Im0) + \
    (Um0 - Im0*Rm)/(KV0*Um0)*N_maxV
Im_maxV = M_maxV*KV0*Um0/(9.55*(Um0 - Im0*Rm)) + Im0

sigma_maxV = (Um_maxV + Im_maxV*Re) / Ubat  # [0, 1]
Ie_maxV = sigma_maxV*Im_maxV  # Missing flag
if Ie_maxV > Ie_max:
    print("ESC current limit has been exceeded!")

Ibat_maxV = n_r*Ie_maxV + Ictrl
if Ibat_maxV > Cbat*Kbat:
    print("Battery current limit has been exceeded!")
Ue_maxV = Ubat - Ibat_maxV*Rbat

t_bat_maxV = (Cbat - Cmin)/Ibat_maxV * 60/1000  # Missing flag
d_maxV = 60*V_maxV*t_bat_maxV

table_maxV = [["Endurance", round(t_bat_maxV, 3), "min"],
              ["Distance", round(d_maxV/1000, 3),
              "km"], ["Fligh Speed", round(V_maxV, 3), "m/s"]]

# print("\nForward Flight:\n")
# print(tabulate(table_maxV, headers=[
#       "Parameter", "Value", "Unit"], tablefmt="github"))

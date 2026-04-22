import numpy as np
from scipy.optimize import fsolve

R = 8.314
T = 520           
dG1 = -40000 
dG2 = -25000 
P_pa = 200000
P_bar = P_pa / 100000

nC2H4_0 = 0.30
nCO0    = 0.30
nH20    = 0.30
nN20    = 0.10

Kp1 = np.exp(-dG1 / (R * T))
Kp2 = np.exp(-dG2 / (R * T))

print(f"Kp1 = {Kp1:.2f}")   
print(f"Kp2 = {Kp2:.2f}")   

def equations(vars):
    x, y = vars

    n_total = 1 - 2*x - y
    
    Kn1 = (x - y) / ((nC2H4_0 - x) * (nCO0 - x) * (nH20 - x - y))
    Kp1_calc = Kn1 * (n_total / P_bar)**2
    
    if x - y <= 0 or (nH20 - x - y) <= 0:
        Kp2_calc = 0
    else:
        Kn2 = y / ((x - y) * (nH20 - x - y))
        Kp2_calc = Kn2 * (n_total / P_bar)
    
    return [Kp1_calc - Kp1, Kp2_calc - Kp2]

x_sol, y_sol = fsolve(equations, (0.1, 0.1))

x, y = x_sol, y_sol
nC2H4  = 0.30 - x
nCO    = 0.30 - x
nH2    = 0.30 - x - y
nC2H5CHO   = x - y
nC3H7OH   = y
nN2    = 0.10
n_total = 1 - 2*x - y

print("\n=== Результаты расчёта ===")
print(f"x (прореагировало в 1 рек) = {x:.6f}")
print(f"y (прореагировало во 2 рек) = {y:.6f}")
print(f"Общее число молей: {n_total:.6f}")

print("\n=== Равновесные мольные доли ===")
print(f"C₂H₄     = {nC2H4/n_total:.6f}")
print(f"CO       = {nCO/n_total:.6f}")
print(f"H₂       = {nH2/n_total:.6f}")
print(f"C₂H₅CHO  = {nC2H5CHO/n_total:.6f}")
print(f"C₃H₇OH   = {nC3H7OH/n_total:.6f}")
print(f"N₂       = {nN2/n_total:.6f}")

alpha_C2H4 = x / 0.30
print(f"\n=== Степень превращения C₂H₄ ===")
print(f"α = {alpha_C2H4:.4f} = {alpha_C2H4*100:.2f}%")

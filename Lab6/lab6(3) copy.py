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

# --- 2. Расчет констант равновесия Ka (Kp) ---
# Считаем по формуле Kp = exp(-dG / RT)
Kp1 = np.exp(-dG1 / (R * T))
Kp2 = np.exp(-dG2 / (R * T))

print(f"Kp1 = {Kp1:.2f}")   
print(f"Kp2 = {Kp2:.2f}")   

# --- 3. Определение системы уравнений ---
def equations(vars):
    x, y = vars
    # x - сколько этилена ушло в первую реакцию (продукт C2H5CHO)
    # y - сколько C2H5CHO ушло во вторую реакцию (продукт C3H7OH)
    
    # Общее число молей из таблицы баланса:
    # (0.3-x) + (0.3-x) + (0.3-x-y) + (x-y) + y + 0.1 = 1 - 2x - y
    n_total = 1 - 2*x - y
    
    # Реакция 1: C2H4 + CO + H2 <=> C2H5CHO
    # Число молей: C2H4: 0.3-x, CO: 0.3-x, H2: 0.3-x-y, Альдегид: x-y
    # Kn1 = n(продукт) / (n1 * n2 * n3)
    Kn1 = (x - y) / ((nC2H4_0 - x) * (nCO0 - x) * (nH20 - x - y))
    # Связь Kp и Kn: Kp = Kn * (P / n_total)^delta_nu. Тут delta_nu = 1 - 3 = -2
    # Kp1 = Kn1 * (P / n_total)**-2  =>  Kp1 = Kn1 * (n_total / P)**2
    Kp1_calc = Kn1 * (n_total / P_bar)**2
    
    # Реакция 2: C2H5CHO + H2 <=> C3H7OH
    # Число молей: Альдегид: x-y, H2: 0.3-x-y, Спирт: y
    # delta_nu = 1 - 2 = -1
    if x - y <= 0 or (nH20 - x - y) <= 0: # Защита от отрицательных концентраций
        Kp2_calc = 0
    else:
        Kn2 = y / ((x - y) * (nH20 - x - y))
        # Kp2 = Kn2 * (P / n_total)**-1 => Kp2 = Kn2 * (n_total / P)
        Kp2_calc = Kn2 * (n_total / P_bar)
    
    return [Kp1_calc - Kp1, Kp2_calc - Kp2]

# --- 4. Решение системы ---
# x0, y0 - начальное приближение. 
# Так как водорода мало (0.3), x+y не может быть больше 0.3. Пробуем малые числа.
x_sol, y_sol = fsolve(equations, (0.1, 0.1))

# --- 5. Расчет состава равновесной смеси ---
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
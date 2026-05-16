import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# --- ИСХОДНЫЕ ДАННЫЕ ВАРИАНТА 17 ---
P = 101325       # Па (давление)
T = 1000         # К (температура)
R = 8.314        # Газовая постоянная
t_end = 0.5      # Время подберем, чтобы увидеть изменение (в секундах)

# Начальные мольные доли
x_CH4 = 0.3
x_H2O = 0.6
x_N2  = 0.1

# Расчет начальных концентраций (C = P*x / R*T) в моль/м3
C_total = P / (R * T)
C_CH4_0 = x_CH4 * C_total
C_H2O_0 = x_H2O * C_total
C_CO_0  = 0.0
C_H2_0  = 0.0
C_CO2_0 = 0.0

# Кинетические параметры (A в 1/с или м3/(моль*с), E в кДж/моль)
# Реакция 1: CH4 + H2O <-> CO + 3H2
A1p, E1p = 1e3, 50
A1m, E1m = 1e2, 124

# Реакция 2: CO + H2O <-> CO2 + H2
A2p, E2p = 2e2, 64
A2m, E2m = 5e2, 130

# Функция для расчета констант скорости k = A * exp(-E / RT)
def calc_k(A, E_kj, T):
    return A * np.exp(-E_kj * 1000 / (R * T))

k1p = calc_k(A1p, E1p, T)
k1m = calc_k(A1m, E1m, T)
k2p = calc_k(A2p, E2p, T)
k2m = calc_k(A2m, E2m, T)

# --- СИСТЕМА ДИФФЕРЕНЦИАЛЬНЫХ УРАВНЕНИЙ ---
def reactor_model(y, t):
    # Концентрации в данный момент времени
    CH4, H2O, CO, H2, CO2 = y
    
    # Скорости элементарных стадий
    # r1: CH4 + H2O <-> CO + 3H2
    r1 = k1p * CH4 * H2O - k1m * CO * (H2**3)
    
    # r2: CO + H2O <-> CO2 + H2
    r2 = k2p * CO * H2O - k2m * CO2 * H2
    
    # Уравнения изменения концентраций (по стехиометрии)
    dCH4_dt = -r1
    dH2O_dt = -r1 - r2
    dCO_dt  =  r1 - r2
    dH2_dt  = 3*r1 + r2
    dCO2_dt =  r2
    
    return [dCH4_dt, dH2O_dt, dCO_dt, dH2_dt, dCO2_dt]

# Вектор начальных условий
y0 = [C_CH4_0, C_H2O_0, C_CO_0, C_H2_0, C_CO2_0]

# Временная сетка
t = np.linspace(0, t_end, 2000)

# Решение системы
solution = odeint(reactor_model, y0, t)
CH4_sol, H2O_sol, CO_sol, H2_sol, CO2_sol = solution.T

# --- ВИЗУАЛИЗАЦИЯ ---
plt.figure(figsize=(10, 6))
plt.plot(t, CH4_sol, 'r-', label='CH₄ (Метан)')
plt.plot(t, H2O_sol, 'b-', label='H₂O (Вода)')
plt.plot(t, CO_sol,  'g-', label='CO (Оксид углерода)')
plt.plot(t, H2_sol,  'm-', label='H₂ (Водород)')
plt.plot(t, CO2_sol, 'k--', label='CO₂ (Углекислый газ)')

plt.xlabel('Время t, с', fontsize=12)
plt.ylabel('Концентрация, моль/м³', fontsize=12)
plt.title(f'Кинетические кривые конверсии метана, T={T}K', fontsize=14)
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

# Печать результатов
print(f"Константы скорости:")
print(f"k1+: {k1p:.4e}, k1-: {k1m:.4e}")
print(f"k2+: {k2p:.4e}, k2-: {k2m:.4e}")
print("\nСостав в конце расчета (моль/м3):")
print(f"CH4: {CH4_sol[-1]:.4f}")
print(f"H2O: {H2O_sol[-1]:.4f}")
print(f"CO:  {CO_sol[-1]:.4f}")
print(f"H2:  {H2_sol[-1]:.4f}")
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

k1 = 0.18
k2 = 0.29

CA0 = 0.028
CB0 = 0.0
CC0 = 0.0
CD0 = 0.0

y0 = [CA0, CB0, CC0, CD0]

def kinetics(y, t):
    CA, CB, CC, CD = y
    r1 = k1 * CA
    r2 = k2 * CB * CC

    dCA_dt = -r1
    dCB_dt = r1 - r2
    dCC_dt = r1 - r2
    dCD_dt = r2
    return [dCA_dt, dCB_dt, dCC_dt, dCD_dt]

t = np.linspace(0, 50, 1000)

solution = odeint(kinetics, y0, t)
CA, CB, CC, CD = solution.T

plt.figure(figsize=(12, 7))
plt.plot(t, CA, 'r-', linewidth=2, label='н-C₇H₁₆ (A)')
plt.plot(t, CB, 'b-', linewidth=2, label='C₃H₆ (B)')
plt.plot(t, CC, 'g-', linewidth=2, label='C₄H₁₀ (C)')
plt.plot(t, CD, 'm-', linewidth=2, label='и-C₄H₁₀ (D)')
plt.xlabel('Время t', fontsize=13)
plt.ylabel('Концентрация, кмоль/м³', fontsize=13)
plt.title('Кинетические кривые\nT = 650 K, k₁ = 0.18, k₂ = 0.29', fontsize=15)
plt.legend(fontsize=12)
plt.grid(True, alpha=0.4)
plt.tight_layout()
plt.show()

print("Конечные концентрации:")
print(f"cA = {CA[-1]:.6f} кмоль/м³")
print(f"cB = {CB[-1]:.6f} кмоль/м³")
print(f"cC = {CC[-1]:.6f} кмоль/м³")
print(f"cD = {CD[-1]:.6f} кмоль/м³")
print(f"Сумма = {CA[-1] + CB[-1] + CC[-1] + CD[-1]:.6f} кмоль/м³")

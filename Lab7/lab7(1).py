import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

print("\nИсходная реакция: 2NOCl -> 2NO + Cl2")

T_target = 490
R_const = 8.314
p_init = 16.81

t_val = np.array([0, 160, 320, 760, 1510, 3820, 5600])
p_total = np.array([16.81, 17.03, 17.24, 17.76, 18.51, 20.10, 20.89])

pA_val = 3 * p_init - 2 * p_total
cA_val = pA_val * 1000 / (R_const * T_target)# По уравнению Менделеева-Клайперона     

print(f"{'t, с':>6}  {'p, кПа':>8}  {'pA, кПа':>9}  {'cA, моль/м³':>12}")
print("-" * 50)
for i in range(len(t_val)):
    print(f"{t_val[i]:6.0f}  {p_total[i]:8.2f}  {pA_val[i]:9.2f}  {cA_val[i]:12.4f}")

y_0 = cA_val
y_1 = np.log(cA_val)
y_2 = 1 / cA_val

slope_0, inter_0, r_0, _, _ = stats.linregress(t_val, y_0)
slope_1, inter_1, r_1, _, _ = stats.linregress(t_val, y_1)
slope_2, inter_2, r_2, _, _ = stats.linregress(t_val, y_2)

R2_val0 = r_0**2
R2_val1 = r_1**2
R2_val2 = r_2**2

R2_results = {0: R2_val0, 1: R2_val1, 2: R2_val2}
best_idx = max(R2_results, key=R2_results.get)
names = {0: 'нулевой', 1: 'первый', 2: 'второй'}

print("\n   Коэффициенты детерминации:")
print(f"Нулевой порядок:  R² = {R2_val0:.4f}")
print(f"Первый порядок:   R² = {R2_val1:.4f}")
print(f"Второй порядок:   R² = {R2_val2:.4f}")
print(f"\nВыбран {names[best_idx]} порядок (R² = {R2_results[best_idx]:.4f})")

plt.figure(figsize=(8, 5))
plt.plot(t_val, cA_val, 'o-', color='darkblue', markersize=7, linewidth=1.5)
plt.xlabel('t, с', fontsize=12)
plt.ylabel('cA, моль/м³', fontsize=12)
plt.title('Кинетическая кривая', fontsize=14)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

t_linspace = np.linspace(0, 5700, 100)

fig, axes = plt.subplots(1, 3, figsize=(18, 6))

axes[0].scatter(t_val, y_0, color='green', s=50, zorder=5)
axes[0].plot(t_linspace, inter_0 + slope_0 * t_linspace, 'k-', linewidth=1.5,
             label=f'R² = {R2_val0:.4f}')
axes[0].set_xlabel('t, с')
axes[0].set_ylabel('cA, моль/м³')
axes[0].set_title('Нулевой порядок')
axes[0].legend()
axes[0].grid(True, alpha=0.3)

axes[1].scatter(t_val, y_1, color='green', s=50, zorder=5)
axes[1].plot(t_linspace, inter_1 + slope_1 * t_linspace, 'k-', linewidth=1.5,
             label=f'R² = {R2_val1:.4f}')
axes[1].set_xlabel('t, с')
axes[1].set_ylabel('ln(cA)')
axes[1].set_title('Первый порядок')
axes[1].legend()
axes[1].grid(True, alpha=0.3)

axes[2].scatter(t_val, y_2, color='green', s=50, zorder=5)
axes[2].plot(t_linspace, inter_2 + slope_2 * t_linspace, 'k-', linewidth=1.5,
             label=f'R² = {R2_val2:.4f}')
axes[2].set_xlabel('t, с')
axes[2].set_ylabel('1/cA')
axes[2].set_title('Второй порядок')
axes[2].legend()
axes[2].grid(True, alpha=0.4)

plt.tight_layout(rect=[0, 0, 1, 0.95]) # Добавляем отступ сверху (0.95) под заголовок
plt.suptitle('Определение порядка реакции (линейные графики)', fontsize=15)
plt.show()

cA_start = cA_val[0]
t1_calc = 3600  # Время из твоего условия

if best_idx == 0:
    k_res = abs(slope_0)
    tau_half = cA_start / (2 * k_res)
    cA_t1 = cA_start - k_res * t1_calc
elif best_idx == 1:
    k_res = abs(slope_1)
    tau_half = np.log(2) / k_res
    cA_t1 = cA_start * np.exp(-k_res * t1_calc)
elif best_idx == 2:
    k_res = abs(slope_2)
    tau_half = 1 / (k_res * cA_start)
    cA_t1 = cA_start / (1 + cA_start * k_res * t1_calc)

deg_conv = (cA_start - cA_t1) / cA_start

print("\n   Итоговые расчеты:")
print(f"Константа скорости k = {k_res:.4e}")
print(f"Время полупревращения τ₁/₂ = {tau_half:.2f} с")
print(f"Концентрация NOCl при t₁ = {t1_calc} с: cA = {cA_t1:.4f} моль/м³")
print(f"Степень превращения при t₁ = {t1_calc} с: α = {deg_conv:.4f}")
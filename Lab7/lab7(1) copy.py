import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

print("\nИсходная реакция: 2NO2 -> 2NO + O2")

T = 643               
R = 8.314             
p0 = 4.88             

t = np.array([0, 30, 60, 120, 180, 240, 320])            
p = np.array([4.88, 5.26, 5.54, 5.91, 6.16, 6.33, 6.49]) 

pA = 3 * p0 - 2 * p               
cA = pA * 1000 / (R * T)          

print(f"{'t, с':>6}  {'p, кПа':>8}  {'pA, кПа':>9}  {'cA, моль/м³':>12}")
print("-" * 50)
for i in range(len(t)):
    print(f"{t[i]:6.0f}  {p[i]:8.2f}  {pA[i]:9.2f}  {cA[i]:12.4f}")

y0 = cA               
y1 = np.log(cA)        
y2 = 1 / cA           

slope0, inter0, r0, _, _ = stats.linregress(t, y0)
slope1, inter1, r1, _, _ = stats.linregress(t, y1)
slope2, inter2, r2, _, _ = stats.linregress(t, y2)

R2_0 = r0**2
R2_1 = r1**2
R2_2 = r2**2

R2_dict = {0: R2_0, 1: R2_1, 2: R2_2}
best_order = max(R2_dict, key=R2_dict.get)
order_names = {0: 'нулевой', 1: 'первый', 2: 'второй'}

print("\n   Коэффициенты детерминации:")
print(f"Нулевой порядок:  R² = {R2_0:.4f}")
print(f"Первый порядок:   R² = {R2_1:.4f}")
print(f"Второй порядок:   R² = {R2_2:.4f}")
print(f"\nВыбран {order_names[best_order]} порядок (R² = {R2_dict[best_order]:.4f})")

plt.figure(figsize=(8, 5))
plt.plot(t, cA, 'o-', color='darkred', markersize=8, linewidth=2)
plt.xlabel('t, с', fontsize=12)
plt.ylabel('cA, моль/м³', fontsize=12)
plt.title('Кинетическая кривая', fontsize=14)
plt.grid(True, alpha=0.4)
plt.tight_layout()
plt.show()

t_line = np.linspace(0, 330, 100)

fig, axes = plt.subplots(1, 3, figsize=(18, 6))

axes[0].scatter(t, y0, color='red', s=60, zorder=5)
axes[0].plot(t_line, inter0 + slope0 * t_line, 'b-', linewidth=2,
             label=f'R² = {R2_0:.4f}')
axes[0].set_xlabel('t, с', fontsize=12)
axes[0].set_ylabel('cA, моль/м³', fontsize=12)
axes[0].set_title('Нулевой порядок\ncA = f(t)', fontsize=13)
axes[0].legend(fontsize=10)
axes[0].grid(True, alpha=0.4)

axes[1].scatter(t, y1, color='red', s=60, zorder=5)
axes[1].plot(t_line, inter1 + slope1 * t_line, 'b-', linewidth=2,
             label=f'R² = {R2_1:.4f}')
axes[1].set_xlabel('t, с', fontsize=12)
axes[1].set_ylabel('ln(cA)', fontsize=12)
axes[1].set_title('Первый порядок\nln(cA) = f(t)', fontsize=13)
axes[1].legend(fontsize=10)
axes[1].grid(True, alpha=0.4)

axes[2].scatter(t, y2, color='red', s=60, zorder=5)
axes[2].plot(t_line, inter2 + slope2 * t_line, 'b-', linewidth=2,
             label=f'R² = {R2_2:.4f}')
axes[2].set_xlabel('t, с', fontsize=12)
axes[2].set_ylabel('1/cA, м³/моль', fontsize=12)
axes[2].set_title('Второй порядок\n1/cA = f(t)', fontsize=13)
axes[2].legend(fontsize=10)
axes[2].grid(True, alpha=0.4)

plt.suptitle('Линейные анаморфозы для определения порядка реакции', fontsize=15, y=1.02)
plt.tight_layout()
plt.show()

cA0 = cA[0]
t1 = 220

if best_order == 0:
    k = abs(slope0)
    tau_half = cA0 / (2 * k)
    cA_t1 = cA0 - k * t1
    if cA_t1 < 0:
        cA_t1 = 0.0
elif best_order == 1:
    k = abs(slope1)
    tau_half = np.log(2) / k
    cA_t1 = cA0 * np.exp(-k * t1)
elif best_order == 2:
    k = abs(slope2)
    tau_half = 1 / (k * cA0)
    cA_t1 = cA0 / (1 + cA0 * k * t1)

alpha = (cA0 - cA_t1) / cA0

print("\n   Расчёты:")
print(f"Константа скорости k = {k:.4e}")
print(f"Время полупревращения τ₁/₂ = {tau_half:.2f} с")
print(f"Концентрация NO₂ при t₁ = {t1} с: cA = {cA_t1:.4f} моль/м³")
print(f"Степень превращения при t₁ = {t1} с: α = {alpha:.4f} ({alpha*100:.1f} %)")

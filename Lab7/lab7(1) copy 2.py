import numpy as np  # Импорт библиотеки для работы с массивами и математикой
import matplotlib.pyplot as plt  # Импорт библиотеки для построения графиков
from scipy import stats  # Импорт модуля для выполнения линейной регрессии

print("\nИсходная реакция: 2NOCl -> 2NO + Cl2")  # Вывод названия реакции

T_target = 490  # Задание температуры в Кельвинах (из варианта 17)
R_const = 8.314  # Задание универсальной газовой постоянной
p_init = 16.81  # Задание начального давления реагента p0 в кПа

t_val = np.array([0, 160, 320, 760, 1510, 3820, 5600])  # Массив временных точек эксперимента
p_total = np.array([16.81, 17.03, 17.24, 17.76, 18.51, 20.10, 20.89])  # Массив общего давления в системе

pA_val = 3 * p_init - 2 * p_total  # Расчет парциального давления NOCl по формуле матбаланса
cA_val = pA_val * 1000 / (R_const * T_target) # Расчет концентрации NOCl (моль/м3) по уравнению Менделеева-Клапейрона     

print(f"{'t, с':>6}  {'p, кПа':>8}  {'pA, кПа':>9}  {'cA, моль/м³':>12}")  # Печать заголовка таблицы
print("-" * 50)  # Печать разделительной линии
for i in range(len(t_val)):  # Цикл по всем экспериментальным точкам
    print(f"{t_val[i]:6.0f}  {p_total[i]:8.2f}  {pA_val[i]:9.2f}  {cA_val[i]:12.4f}")  # Печать строки таблицы данных

y_0 = cA_val  # Координаты Y для проверки нулевого порядка (c)
y_1 = np.log(cA_val)  # Координаты Y для проверки первого порядка (ln c)
y_2 = 1 / cA_val  # Координаты Y для проверки второго порядка (1/c)

slope_0, inter_0, r_0, _, _ = stats.linregress(t_val, y_0)  # Регрессия для нулевого порядка
slope_1, inter_1, r_1, _, _ = stats.linregress(t_val, y_1)  # Регрессия для первого порядка
slope_2, inter_2, r_2, _, _ = stats.linregress(t_val, y_2)  # Регрессия для второго порядка

R2_val0 = r_0**2  # Коэффициент детерминации для нулевого порядка
R2_val1 = r_1**2  # Коэффициент детерминации для первого порядка
R2_val2 = r_2**2  # Коэффициент детерминации для второго порядка

R2_results = {0: R2_val0, 1: R2_val1, 2: R2_val2}  # Словарь с результатами R2
best_idx = max(R2_results, key=R2_results.get)  # Поиск индекса самого точного (максимального) R2
names = {0: 'нулевой', 1: 'первый', 2: 'второй'}  # Словарь названий порядков

print("\n   Коэффициенты детерминации:")  # Заголовок вывода R2
print(f"Нулевой порядок:  R² = {R2_val0:.4f}")  # Вывод R2 для 0 порядка
print(f"Первый порядок:   R² = {R2_val1:.4f}")  # Вывод R2 для 1 порядка
print(f"Второй порядок:   R² = {R2_val2:.4f}")  # Вывод R2 для 2 порядка
print(f"\nВыбран {names[best_idx]} порядок (R² = {R2_results[best_idx]:.4f})")  # Вывод итогового выбранного порядка

plt.figure(figsize=(8, 5))  # Создание окна для первого графика
plt.plot(t_val, cA_val, 'o-', color='darkblue', markersize=7, linewidth=1.5)  # Отрисовка экспериментальных точек c(t)
plt.xlabel('t, с', fontsize=12)  # Подпись оси X
plt.ylabel('cA, моль/м³', fontsize=12)  # Подпись оси Y
plt.title('Кинетическая кривая', fontsize=14)  # Заголовок графика
plt.grid(True, alpha=0.3)  # Включение сетки на графике
plt.tight_layout()  # Авто-выравнивание элементов графика
plt.show()  # Отображение графика на экране

t_linspace = np.linspace(0, 5700, 100)  # Создание 100 ровных точек времени для отрисовки линий регрессии
fig, axes = plt.subplots(1, 3, figsize=(18, 6))  # Создание окна с 3-мя подграфиками

axes[0].scatter(t_val, y_0, color='green', s=50, zorder=5)  # Точки для нулевого порядка
axes[0].plot(t_linspace, inter_0 + slope_0 * t_linspace, 'k-', linewidth=1.5, label=f'R² = {R2_val0:.4f}')  # Линия для 0 порядка
axes[0].set_xlabel('t, с'); axes[0].set_ylabel('cA, моль/м³'); axes[0].set_title('Нулевой порядок'); axes[0].legend(); axes[0].grid(True, alpha=0.3)  # Настройка осей

axes[1].scatter(t_val, y_1, color='green', s=50, zorder=5)  # Точки для первого порядка
axes[1].plot(t_linspace, inter_1 + slope_1 * t_linspace, 'k-', linewidth=1.5, label=f'R² = {R2_val1:.4f}')  # Линия для 1 порядка
axes[1].set_xlabel('t, с'); axes[1].set_ylabel('ln(cA)'); axes[1].set_title('Первый порядок'); axes[1].legend(); axes[1].grid(True, alpha=0.3)  # Настройка осей

axes[2].scatter(t_val, y_2, color='green', s=50, zorder=5)  # Точки для второго порядка
axes[2].plot(t_linspace, inter_2 + slope_2 * t_linspace, 'k-', linewidth=1.5, label=f'R² = {R2_val2:.4f}')  # Линия для 2 порядка
axes[2].set_xlabel('t, с'); axes[2].set_ylabel('1/cA'); axes[2].set_title('Второй порядок'); axes[2].legend(); axes[2].grid(True, alpha=0.4)  # Настройка осей

plt.tight_layout(rect=[0, 0, 1, 0.95])  # Выравнивание графиков с местом под общий заголовок
plt.suptitle('Определение порядка реакции (линейные графики)', fontsize=15)  # Общий заголовок окна
plt.show()  # Отображение графиков анаморфоз

cA_start = cA_val[0]  # Начальная концентрация в момент времени t=0
t1_calc = 3600  # Время t1 для расчета параметров из условия задачи

if best_idx == 0:  # Если выбран нулевой порядок
    k_res = abs(slope_0)  # Константа скорости k
    tau_half = cA_start / (2 * k_res)  # Время полупревращения
    cA_t1 = cA_start - k_res * t1_calc  # Концентрация в момент t1
elif best_idx == 1:  # Если выбран первый порядок
    k_res = abs(slope_1)  # Константа скорости k
    tau_half = np.log(2) / k_res  # Время полупревращения
    cA_t1 = cA_start * np.exp(-k_res * t1_calc)  # Концентрация в момент t1
elif best_idx == 2:  # Если выбран второй порядок
    k_res = abs(slope_2)  # Константа скорости k
    tau_half = 1 / (k_res * cA_start)  # Время полупревращения
    cA_t1 = cA_start / (1 + cA_start * k_res * t1_calc)  # Концентрация в момент t1

deg_conv = (cA_start - cA_t1) / cA_start  # Расчет степени превращения (конверсии)

p_final = p_total[-1]  # Конечное общее давление
x_final = p_final - p_init  # Давление превратившейся части хлора
pA_final = pA_val[-1]  # Конечное давление NOCl
pNO_final = 2 * x_final  # Конечное давление NO
pCl2_final = x_final  # Конечное давление Cl2

cNOCl_eq = (pA_final * 1000) / (R_const * T_target)  # Равновесная концентрация NOCl
cNO_eq = (pNO_final * 1000) / (R_const * T_target)  # Равновесная концентрация NO
cCl2_eq = (pCl2_final * 1000) / (R_const * T_target)  # Равновесная концентрация Cl2

Kc = (cNO_eq**2 * cCl2_eq) / (cNOCl_eq**2)  # Расчет константы равновесия Kc

if best_idx == 0:  # Определение единиц измерения k
    k_unit = "моль / (м³ * с)"
elif best_idx == 1:
    k_unit = "с⁻¹"
else:
    k_unit = "м³ / (моль * с)"

print("\n   Итоговые расчеты:")  # Вывод финальных результатов в консоль
print(f"Константа скорости k = {k_res:.4e} {k_unit}") # Печать k
print(f"Время полупревращения τ₁/₂ = {tau_half:.2f} с") # Печать времени полураспада
print(f"Концентрация NOCl при t₁ = {t1_calc} с: cA = {cA_t1:.4f} моль/м³") # Печать c(t1)
print(f"Конверсия при t₁ = {t1_calc} с: α = {deg_conv:.4f} ({deg_conv*100:.2f}%)") # Печать конверсии
print(f"Константа равновесия Kc = {Kc:.4f}") # Печать константы равновесия
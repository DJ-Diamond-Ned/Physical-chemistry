import numpy as np
from scipy.optimize import fsolve

T = 1000
R = 8.314462618

coef_H2 = [0.23443029E+01, 0.79804248E-02, -0.19477917E-04, 0.20156967E-07, -0.73760289E-11, -0.91792413E+03, 0.68300218E+00]
coef_CO = [0.35795335E+01, -0.61035369E-03, 0.10168143E-05, 0.90700586E-09, -0.90442449E-12, -0.14344086E+05, 0.35084093E+01]
coef_H2O = [7.25575005E+01, -6.62445402E-01, 2.56198746E-03, -4.36591923E-06, 2.78178981E-09, -4.18865499E+04, -2.88280137E+02]
coef_CO2 = [0.23568130E+01, 0.89841299E-02, -0.71220632E-05, 0.24573008E-08, -0.14288548E-12, -0.48371971E+05, 0.99009035E+01]

def calc_H(coef, T):
    return R * T * (coef[0] + coef[1]*T/2 + coef[2]*pow(T,2)/3 + coef[3]*pow(T,3)/4 + coef[4]*pow(T,4)/5 + coef[5]/T)

def calc_S(coef, T):
    return (coef[0] * np.log(T) + coef[1]*T + (coef[2]/2)*T**2 + (coef[3]/3)*T**3 + (coef[4]/4)*T**4 + coef[6]) * R
     

print("1. ΔG⁰(Т) по полиномам NASA")

H_H2  = calc_H(coef_H2,  T)
S_H2  = calc_S(coef_H2,  T)
H_CO = calc_H(coef_CO, T)
S_CO= calc_S(coef_CO, T)
H_H2O  = calc_H(coef_H2O,  T)
S_H2O  = calc_S(coef_H2O,  T)
H_CO2  = calc_H(coef_CO2,  T)
S_CO2  = calc_S(coef_CO2,  T)

delta_H = -1 * H_CO - 1 * H_H2O +  1 * H_CO2 + 1 * H_H2
delta_S = -1 * S_CO - 1 * S_H2O +  1 * S_CO2 + 1 * S_H2
delta_G = delta_H - T * delta_S

print(f"ΔH⁰ = {delta_H:.2f} Дж/моль = {delta_H/1000:.2f} кДж/моль")
print(f"ΔS⁰ = {delta_S:.2f} Дж/(моль·К)")
print(f"ΔG⁰ = {delta_G:.2f} Дж/моль = {delta_G/1000:.2f} кДж/моль")

print("\n2. Мольные доли веществ в начале реакции")

P_H2_0_Pa = 1.6 * pow(10,4)
P_CO_0_Pa  = 5 * pow(10,4)
P_H2O_0_Pa  = 2 * pow(10,4)
P_CO2_0_Pa  = 0.7 * pow(10,4)

P_H2_0 = P_H2_0_Pa / 100000
P_CO_0  = P_CO_0_Pa / 100000
P_H2O_0  = P_H2O_0_Pa / 100000
P_CO2_0  = P_CO2_0_Pa / 100000

P_total_0 = P_H2_0 + P_CO_0 + P_H2O_0 + P_CO2_0

y_H2_0 = P_H2_0 / P_total_0
y_CO_0  = P_CO_0  / P_total_0
y_H2O_0  = P_H2O_0  / P_total_0
y_CO2_0  = P_CO2_0  / P_total_0

print(f"Общее давление в начале: P_total⁰ = {P_total_0:.4f} бар")
print(f"Мольные доли в начале:")
print(f"  y_H2⁰ = {y_H2_0:.4f}")
print(f"  y_CO⁰  = {y_CO_0:.4f}")
print(f"  y_H2O⁰  = {y_H2O_0:.4f}")
print(f"  y_CO2⁰  = {y_CO2_0:.4f}")
print(f"Сумма = {y_H2_0 + y_CO_0 + y_H2O_0 + y_CO2_0:.4f}")

print("\n3. Константа равновесия Keq ")

Keq = np.exp(-delta_G / (R * T))
print(f"Keq = {Keq:.4e}")

print("\n4. Равновесные концентрации всех веществ в системе")

print("\n Таблица (в мольных долях, Σn₀ = 1):")
print("=" * 73)
print(f"| Компонент |     H₂      |     CO      |     H₂O      |     CO₂        |")
print("=" * 73)
print(f"| Было      |   {y_H2_0:.4f}    |   {y_CO_0:.4f}    |   {y_H2O_0:.4f}     |     {y_CO2_0:.4f}     |")
print(f"| Изменение |     +x      |     -x      |    -x        |    +x          |")
print(f"| Стало     | {y_H2_0:.4f} + x  | {y_CO_0:.4f} - x  | {y_H2O_0:.4f} - x   |     {y_CO2_0:.4f} + x |")
print("=" * 73)
print(f"Σn = 1")
print(f"Pобщ = Σn * Pобщ⁰ = 1  * {P_total_0:.1f} бар")
print(f"Kp = (CO2 * H2) / (CO * H2O)")

def equation(x):
    if x >= y_H2O_0 or x >= y_CO_0:
        return 1e10

    p_H2  = P_H2_0 + x
    p_CO  = P_CO_0 - x
    p_H2O = P_H2O_0 - x
    p_CO2 = P_CO2_0 + x

    Kp_expr = (p_CO2 * p_H2) / (p_CO * p_H2O)
    return Kp_expr - Keq

x_sol = fsolve(equation, 1e-6)[0]

print(f"\nСтепень превращения x = {x_sol:.6f}")

p_H2_eq  = P_H2_0 + x_sol
p_CO_eq  = P_CO_0 - x_sol
p_H2O_eq = P_H2O_0 - x_sol
p_CO2_eq = P_CO2_0 + x_sol

print("Равновесные концентрации")
print(f"p(H2)  = {p_H2_eq:.4f} бар")
print(f"p(CO)  = {p_CO_eq:.4f} бар")
print(f"p(H2O) = {p_H2O_eq:.4f} бар")
print(f"p(CO2) = {p_CO2_eq:.4f} бар")

print("\n5. Равновесная степень превращения CO")
alpha = x_sol / y_CO_0
print(f"α = {alpha:.6f} = {alpha*100:.4f}%")

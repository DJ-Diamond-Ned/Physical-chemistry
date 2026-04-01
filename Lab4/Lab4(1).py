#17.	Этилацетат + вода
# 
import math
from scipy.optimize import minimize

R = 8.314
T = 298.15
Bar = 750.062

A1 = 16.15
B1 = 2790.50
C1 = -577.15
A2 = 18.30
B2 = 3816.44
C2 = -46.13

X1 = [0.1, 0.3, 0.5, 0.7, 0.9]
X2 = [1 - x for x in X1]
Y1 = [0.796, 0.770, 0.813, 0.774, 0.839]
Y2 = [1 - y for y in Y1]
P  = [0.1397, 0.1427, 0.1438, 0.1451, 0.1444]

pBar1 = 38.30
pBar2 = 220.48

Pv1 = (math.exp((A1-B1)/(T + C1)))/Bar
Pv2 = (math.exp((A2-B2)/(T + C2)))/Bar

g1 = [Y1[i] * P[i] / (X1[i] * Pv1) for i in range(5)]
g2 = [Y2[i] * P[i] / (X2[i] * Pv2) for i in range(5)]

gE_exp = [R * T * (X1[i] * math.log(g1[i]) + X2[i] * math.log(g2[i])) for i in range(5)]

def gE_wilson(x1, L12, L21):
    x2 = 1 - x1
    return -R * T * (x1 * math.log(x1 + L12 * x2) + x2 * math.log(x2 + L21 * x1))

def objective(params):
    L12, L21 = params
    if L12 <= 0 or L21 <= 0:
        return 1e10
    return sum(abs(gE_exp[i] - gE_wilson(X1[i], L12, L21)) for i in range(5))

def main():
    result = minimize(objective, [1.0,1.0], method="Nelder-Mead", options={"xatol": 1e-10, "fatol": 1e-10, "maxiter": 2000})

    L12, L21 = result.x
    print(f"Λ12 = {L12:.6f}")
    print(f"Λ21 = {L21:.6f}")

if __name__ == "__main__":
    main()

import math
#17.	Этилацетат + вода
R = 8.314462618
A1 = 16.15
B1 = 2790.50
C1 = -577.15
A2 = 18.30
B2 = 3816.44
C2 = -46.13
T = 298.15
Bar = 750.062
pBar1 = 38.30
pBar2 = 220.48

P10 = (math.exp((A1-B1)/(T + C1)))/Bar
P20 = (math.exp((A2-B2)/(T + C2)))/Bar

x1 = [0.1, 0.3, 0.5, 0.7, 0.9]
y1 = [0.796, 0.770, 0.813, 0.774, 0.839]
x2 = []
y2 = []

for i in range(5):
    x2.append(1 - x1[i])
    y2.append(1 - y1[i])

gamma1 = []
gamma2 = []
ge = []

def calc_Gamma(y, x, pBar, P):
    return (y*pBar)/(x*pBar*P)

def calc_Ge(x1, x2, gamma1, gamma2):
    return R * T * (x1 * math.log(gamma1) + x2 * math.log(gamma2))

for x, y in zip(x1, y1):
    gamma1.append(calc_Gamma(y, x, pBar1, P10))

for x, y in zip(x2, y2):
    gamma2.append(calc_Gamma(y, x, pBar2, P20))
    
for x1, x2, gamma1, gamma2 in zip(x1, x2, gamma1, gamma2):
    ge.append(calc_Ge(x1,x2,gamma1,gamma2))

for i in range(len(ge)):
    print(f"Значение {i + 1} точки = {ge[i]}")








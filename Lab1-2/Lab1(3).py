#O2 Ar
import math
T = 298
P = 1.01 * 10**5
V1 = 1 * 10 ** -4
V2 = 3 * 10 ** -4
R = 8.31
n1 = (P*V1)/(R*T)
n2 = (P*V2)/(R*T)
V = V1 + V2
H = ((n1 * R * math.log(V/V1)) + (n2 * R * math.log(V/V2)))
print(f"Изменение энтропии при смешении = {H:.3f} Дж/к")

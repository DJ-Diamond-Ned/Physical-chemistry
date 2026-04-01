import math

T = 340
R = 8.314462618

coef_ClC3H4 = [3.16995308E+00, 1.57436972E-02, 1.85511623E-05, -3.88261489E-08, 1.77294125E-11, 1.47902471E+04, 1.32175604E+01]
coef_O2 = [3.78245636E+00, -2.99673415E-03, 9.84730200E-06, -9.68129508E-09, 3.24372836E-12, -1.06394356E+03, 3.65767573E+00]
coef_CO2 = [0.23568130E+01, 0.89841299E-02, -0.71220632E-05, 0.24573008E-08, -0.14288548E-12, -0.48371971E+05, 0.99009035E+01]
coef_H2O = [0.41986352E+01, -0.20364017E-02, 0.65203416E-05, -0.54879269E-08, 0.17719680E-11, -0.30293726E+05, -0.84900901E+00]
coef_HCl = [0.34637647E+01, 0.47648423E-03, -0.20030122E-05, 0.33171437E-08, -0.14495818E-11, -0.12144352E+05, 0.26642828E+01]

def calc_H(coef, T):
    return R * T * (coef[0] + coef[1]*T/2 + coef[2]*pow(T,2)/3 + coef[3]*pow(T,3)/4 + coef[4]*pow(T,4)/5 + coef[5]/T)

def calc_S(coef, T):
    return R * (coef[0] * math.log(T) + coef[1] * T + (coef[2]/2) * pow(T,2) + (coef[3]/3) * pow(T,3) + (coef[4]/4) * pow(T,4) + coef[6])

def calc_Cp(coef, T):
    return R * (coef[0] + coef[1]*T + coef[2]*pow(T,2) + coef[3]*pow(T,3) + coef[4]*pow(T,4))

H = calc_H(coef_ClC3H4, T)
S = calc_S(coef_ClC3H4, T)
Cp = calc_Cp(coef_ClC3H4, T)

H_ClC3H4 = calc_H(coef_ClC3H4, T)
H_O2 = calc_H(coef_O2, T)
H_CO2 = calc_H(coef_CO2, T)
H_H2O = calc_H(coef_H2O, T)
H_HCl = calc_H(coef_HCl, T)

S_ClC3H4 = calc_S(coef_ClC3H4, T)
S_O2 = calc_S(coef_O2, T)
S_CO2 = calc_S(coef_CO2, T)
S_H2O = calc_S(coef_H2O, T)
S_HCl = calc_S(coef_HCl, T)

delta_H = 6*H_CO2 + 2*H_HCl + 2*H_H2O - (2*H_ClC3H4 + 7*H_O2)
delta_S = 6*S_CO2 + 2*S_HCl + 2*H_H2O - (2*S_ClC3H4 + 7*S_H2O)
delta_G = delta_H - T * delta_S

print(f"Тепловой эффект реакции:   {delta_H:.2f} Дж/моль")
print(f"Энтропия ClC3H4 : {S:>.2f} Дж/(моль*К)")
print(f"Изменение энергии Гиббса: {delta_G:>.2f} Дж/моль")
print(f"Теплоемкость ClC3H4 : {Cp:>.2f} Дж/моль")

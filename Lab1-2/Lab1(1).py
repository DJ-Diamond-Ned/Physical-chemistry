import math

#Коэфициенты полинома
a1 = 1.25245480E+01
a2 = -1.01018826E-02
a3 = 2.21992610E-04
a4 = -2.84863722E-07
a5 = 1.12410138E-10

#Газовая постоянная
R = 8.314462618

#Расчет полинома
def calc(T):
    return R * (a1 + a2*T + a3*pow(T,2) + a4*pow(T,3) + a5*pow(T,4))

#Данные из сайта(константы)
t_NIST = [298.15, 300, 400, 500, 600]
c_NIST = [187.8, 188.7, 239.74, 286.81, 326.77]

print("T(K)\tcp_Расчетное\tcp_NIST")

#Подстановка данных из файла
for T, cp_n in zip(t_NIST, c_NIST):
    cp_c = calc(T)
    print(f"{T}\t{cp_c:.2f}\t\t{cp_n}")

T_for_err = 400
c_err_calc_400 = calc(T_for_err)
c_err_NIST_400 = 286.81

absolute_error = abs(c_err_calc_400 - c_err_NIST_400)
relative_error = (absolute_error/c_err_NIST_400)*100


print(f"Абсолютная ошибка = {absolute_error:.2f} Дж/(кг*Моль)")
print(f"Относительная ошибка при 400К = {relative_error:.2f}%")

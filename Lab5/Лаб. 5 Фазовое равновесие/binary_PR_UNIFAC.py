from thermo import *
from thermo.unifac import DOUFSG, DOUFIP2016
from dataclasses import dataclass
# Load constants and properties
constants, properties = ChemicalConstantsPackage.from_IDs(["acetone", "chloroform"])
# Objects are initialized at a particular condition

@dataclass(slots=True, frozen=True)
class CONSTANTS:
    T = 298.15
    P = 1e5
    zs = [0.4, 0.6]

# Use Peng-Robinson for the vapor phase
eos_kwargs = {"Pcs": constants.Pcs, "Tcs": constants.Tcs, "omegas": constants.omegas}
gas = CEOSGas(PRMIX, HeatCapacityGases = properties.HeatCapacityGases, eos_kwargs = eos_kwargs)

# Configure the activity model
GE = UNIFAC.from_subgroups(chemgroups = constants.UNIFAC_Dortmund_groups, version = 1, T = CONSTANTS.T, xs = CONSTANTS.zs, interaction_data = DOUFIP2016, subgroups = DOUFSG)
# Configure the liquid model with activity coefficients
liquid = GibbsExcessLiquid(VaporPressures = properties.VaporPressures, HeatCapacityGases = properties.HeatCapacityGases,
    VolumeLiquids = properties.VolumeLiquids, GibbsExcessModel = GE, equilibrium_basis = "Psat", caloric_basis = "Psat", T = CONSTANTS.T, P = CONSTANTS.P, zs = CONSTANTS.zs)

# Create a flasher instance, assuming only vapor-liquid behavior
flasher = FlashVL(constants, properties, liquid = liquid, gas = gas)

# Create a T-xy plot at P bar
_ = flasher.plot_Txy(P = CONSTANTS.P, pts = 100)

# Create a P-xy plot at T Kelvin
_ = flasher.plot_Pxy(T = CONSTANTS.T, pts = 100)

# Create a xy diagram at T Kelvin
_ = flasher.plot_xy(T = CONSTANTS.T, pts = 100)

liquid2 = CEOSLiquid(PRMIX, eos_kwargs, HeatCapacityGases = properties.HeatCapacityGases, T = CONSTANTS.T, P = CONSTANTS.P, zs = CONSTANTS.zs)
flasher2 = FlashVLN(constants, properties, liquids = [liquid, liquid2], gas = gas)
res = flasher2.flash(T = CONSTANTS.T, P = CONSTANTS.P, zs = CONSTANTS.zs)
print(f"При температуре {CONSTANTS.T} K и давлении {CONSTANTS.P / 1e5} бар присутствуют ровно {res.phase_count} фазы")

if (res.VF > 0):
    print(res.gas.zs)
     
if (res.VF == 1):
    print("Только газ") 
     
else:
    print(f"Liquid0: {res.liquid0.zs}") 
    if (res.liquid_count > 1):
        print("LIQUID PHASE SEPARATION")
            


@dataclass(slots=True, frozen=True)
class CONSTANTS_FOR_ERRORS:
    T = 298.15
    Pressure_bar = 1.01

    experimental_points_from_vlecalc = [
        {'T': 333.640, 'x_exp': 0.73,     'y_exp': 0.817000},
        {'T': 332.815, 'x_exp': 0.78,     'y_exp': 0.861500},
        {'T': 333.320, 'x_exp': 0.75,     'y_exp': 0.835822},
        {'T': 332.480, 'x_exp': 0.80,     'y_exp': 0.878659},
        {'T': 331.614, 'x_exp': 0.85,     'y_exp': 0.916240}
    ]

Pressure_pascal = CONSTANTS_FOR_ERRORS.Pressure_bar * 1e5
errors_x, errors_y = [], []

print("\t\t\t\t\t\t\t\t_____ ----- РАСЧЁТ АБСОЛЮТНЫХ ОШИБОК ----- _____")
print("\n\t\tN \t  T(K) \t\tx_exp \t\t y_exp \t\tz=(x+y)/2 \tx_pred \t\t y_pred\t\t |Δx| \t\t  |Δy|")

print("-" * (69 + 67 + 52))

for i, point in enumerate(CONSTANTS_FOR_ERRORS.experimental_points_from_vlecalc, 1):
    z_calc = (point["x_exp"] + point["y_exp"]) / 2
    zs_point = [z_calc, 1 - z_calc]
    
    res_point = flasher.flash(T = point["T"], P = Pressure_pascal, zs = zs_point)
    
    if (res_point.liquid_count > 0):
        x_pred = res_point.liquid0.zs[0]
    else:
        x_pred = 0.0
    
    if ( (res_point.VF > 0) and (res_point.gas is not None) ):
        y_pred = res_point.gas.zs[0]
    else:
        y_pred = 0.0
    
    delta_x = abs(x_pred - point["x_exp"])
    delta_y = abs(y_pred - point["y_exp"])
    
    errors_x.append(delta_x)
    errors_y.append(delta_y)
    
    phase_info = f"({res_point.phase_count} фазы, VF = {res_point.VF:.2f})"    
    print(f"\t\t{i} \t{point["T"]:.4f} \t{point["x_exp"]:.4f} \t\t{point["y_exp"]:.6f} \t{z_calc:.6f} \t{x_pred:.4f} \t\t{y_pred:.6f} \t{delta_x:.5f} \t{delta_y:.7f} \t{phase_info}")

MAE_x = sum(errors_x) / len(errors_x)
MAE_y = sum(errors_y) / len(errors_y)

print(f"\nСредняя абсолютная ошибка (MAE) для жидкой фазы (x): {MAE_x:.6f}")
print(f"Средняя абсолютная ошибка (MAE) для паровой фазы (y): {MAE_y:.6f}")

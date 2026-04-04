from thermo import ChemicalConstantsPackage, CEOSGas, CEOSLiquid, PRMIX, FlashVL, FlashVLN
from thermo.interaction_parameters import IPDB
# Load constants and properties
constants, properties = ChemicalConstantsPackage.from_IDs(['acetone', 'chloroform'])
# Objects are initialized at a particular condition
T = 298.15
P = 1e5
zs = [.5, .5]

# Use Peng-Robinson for both the vapor and the liquid phases
k12 = 0.0159


kijs = [[0, k12],
        [k12, 0]]
print(k12)
eos_kwargs = dict(Tcs=constants.Tcs, Pcs=constants.Pcs, omegas=constants.omegas, kijs=kijs)
gas = CEOSGas(PRMIX, eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases, T=T, P=P, zs=zs)
liquid = CEOSLiquid(PRMIX, eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases, T=T, P=P, zs=zs)
#gas = CEOSGas(PRMIX, eos_kwargs=eos_kwargs)
#liquid = CEOSLiquid(PRMIX, eos_kwargs=eos_kwargs)

# Create a flasher instance, assuming only vapor-liquid behavior
flasher = FlashVL(constants, properties, liquid=liquid, gas=gas)

# Create a T-xy plot at P bar
_ = flasher.plot_Txy(P=P, pts=100)

# Create a P-xy plot at T Kelvin
_ = flasher.plot_Pxy(T=T, pts=100)

# Create a xy diagram at T Kelvin
_ = flasher.plot_xy(T=T, pts=100)

liquid2 = CEOSLiquid(PRMIX, eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases, T=T, P=P, zs=zs)
flasher2 = FlashVLN(constants, properties, liquids=[liquid, liquid2], gas=gas)
res = flasher2.flash(T=T, P=P, zs=zs)
print('There are %s phases present at %f K and %f bar' %(res.phase_count,T,P/1e5))
if res.VF > 0:
	print(res.gas.zs)
if res.VF == 1:  # Есть только пар
	print("Only vapour")
else:
	print("Liquid0: ")
	print(res.liquid0.zs)
	if res.liquid_count>1:
            print("LIQUID PHASE SEPARATION")
		#print("Liquid1: ")
		#print(res.liquid1.zs)

print("Рассчет абсолютных ошибок")

experimental_points = [
    {'T': 334.817, 'x_exp': 0.06,     'y_exp': 0.037680},
    {'T': 336.962, 'x_exp': 0.30,     'y_exp': 0.282070},
    {'T': 336.061, 'x_exp': 0.55,     'y_exp': 0.614500},
    {'T': 333.320, 'x_exp': 0.75,     'y_exp': 0.835822},
    {'T': 330.390, 'x_exp': 0.92,     'y_exp': 0.960424}
]

P_bar = 1.01  
P_Pa = P_bar * 1e5 

errors_x = []
errors_y = [] 

print("\n{:<4} {:<10} {:<10} {:<10} {:<10} {:<12} {:<12} {:<12} {:<12}".format(
    "N", "T(K)", "x_exp", "y_exp", "z=(x+y)/2", "x_pred", "y_pred", "|Δx|", "|Δy|"))
print("-" * 95)

for i, point in enumerate(experimental_points, 1):
    z_calc = (point['x_exp'] + point['y_exp']) / 2
    zs_point = [z_calc, 1 - z_calc]  # [ацетонитрил, метанол]
    
    res_point = flasher2.flash(T=point['T'], P=P_Pa, zs=zs_point)
    
    if res_point.liquid_count > 0:
        x_pred = res_point.liquid0.zs[0]  
    else:
        x_pred = 0.0 
    
    if res_point.VF > 0 and res_point.gas is not None:
        y_pred = res_point.gas.zs[0] 
    else:
        y_pred = 0.0  
    
    delta_x = abs(x_pred - point['x_exp'])
    delta_y = abs(y_pred - point['y_exp'])
    
    errors_x.append(delta_x)
    errors_y.append(delta_y)
    
    print("{:<4} {:<10.3f} {:<10.4f} {:<10.6f} {:<10.4f} {:<12.4f} {:<12.6f} {:<12.6f} {:<12.6f}".format(
        i, point['T'], point['x_exp'], point['y_exp'], z_calc, x_pred, y_pred, delta_x, delta_y))

MAE_x = sum(errors_x) / len(errors_x)
MAE_y = sum(errors_y) / len(errors_y)

print("\тСредняя абсолютная ошибка (MAE) для жидкой фазы (x): {:.6f}".format(MAE_x))
print("Средняя абсолютная ошибка (MAE) для паровой фазы (y): {:.6f}".format(MAE_y))
print("="*100)

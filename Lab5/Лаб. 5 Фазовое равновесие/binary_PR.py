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
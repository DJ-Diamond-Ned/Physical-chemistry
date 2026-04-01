from thermo import *
from thermo.unifac import DOUFSG, DOUFIP2016
# Load constants and properties
constants, properties = ChemicalConstantsPackage.from_IDs(['acetone', 'chloroform'])
# Objects are initialized at a particular condition
T = 298.15
P = 1e5
zs = [.4, .6]

# Use Peng-Robinson for the vapor phase
eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas}
gas = CEOSGas(PRMIX, HeatCapacityGases=properties.HeatCapacityGases, eos_kwargs=eos_kwargs)

# Configure the activity model
GE = UNIFAC.from_subgroups(chemgroups=constants.UNIFAC_Dortmund_groups, version=1, T=T, xs=zs,
                           interaction_data=DOUFIP2016, subgroups=DOUFSG)
# Configure the liquid model with activity coefficients
liquid = GibbsExcessLiquid(
    VaporPressures=properties.VaporPressures,
    HeatCapacityGases=properties.HeatCapacityGases,
    VolumeLiquids=properties.VolumeLiquids,
    GibbsExcessModel=GE,
    equilibrium_basis='Psat', caloric_basis='Psat',
    T=T, P=P, zs=zs)

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
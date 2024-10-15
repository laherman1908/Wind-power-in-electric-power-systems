import numpy as np
from matplotlib import pyplot as plt

# Task 2 - Surge impedance loading and max transmitted power for a three-phase cable

# # Leser inn parametre fra task 2
Vn = 150e3          # Nominal voltage [V]
In = 1088           # Nominal current [A]
Rac = 0.0205        # Ohm/km
L = 0.352e-3        # Line inductance [H/km]
C = 0.233e-6        # Line capacitance
fn = 50             # Nominal frequency

w = 2*np.pi*fn
B = w*C
X = w*L


# # Leser inn parametre fra paper (https://sintef.brage.unit.no/sintef-xmlui/bitstream/handle/11250/2429845/Vrana2016oov.pdf?sequence=2&isAllowed=y)
# Vn = 150e3          # Nominal voltage [V]
# In = 825           # Nominal current [A]
# Rac = 0.0275        # Ohm/km
# L = 0.39e-3        # Line inductance [H/km]
# C = 0.18e-6        # Line capacitance [F/km]
# fn = 50             # Nominal frequency

# w = 2*np.pi*fn
# B = w*C
# X = w*L

# # Leser inn parametre fra Offshore HVAC Electric Power Transmission Systems – an Example
# Vn = 150e3          # Nominal voltage [V]
# In = 1088         # Nominal current [A]
# Rac = 0.0205        # Ohm/km
# L = 0.352e-3        # Line inductance [H/km]
# C = 0.233e-6        # Line capacitance [F/km]
# fn = 50             # Nominal frequency

# w = 2*np.pi*fn
# B = w*C
# X = w*L

def max_power_transfer(Vn : float) -> dict:
    
    Pmax = np.sqrt(3) * Vn * In
    lmax = np.sqrt(12)*In/(B*Vn)

    print(f"lmax: {lmax}")

    Kloss = (3*Rac*lmax*In**2) / (Pmax)

    res = {"power" : [],
           "length" : [],
           "Kloss" : Kloss}
    
    res["Kloss"] = Kloss

    for l in range(1, int(round(lmax))+1):
        ls = l/lmax
        sqrt_term = 1 - (ls**2)
        if sqrt_term < 0:
            sqrt_term = 0
        Ptrans = Pmax*(np.sqrt(sqrt_term)- Kloss*(ls - (2/3)*ls**3))

        res["power"].append(Ptrans*1e-6)        # Convert to MW
        res["length"].append(l)

    return res


def max_power_transfer_2(Vn : float) -> dict[str, list[float]]:

    Smax = np.sqrt(3)*Vn*In
    res = {"power" : [],
        "length" : []
        }

    for l in range(0,400):
        sqrt_term = Smax**2 - ((3*B*l)*(Vn/np.sqrt(3))**2)**2
        if sqrt_term < 0:
            sqrt_term = 0
        print(f"sqrt_term: {sqrt_term}")
        Pmax  = np.sqrt(sqrt_term)
        res["power"].append(Pmax*1e-6)
        res["length"].append(l)
        if Pmax <= 0:
            break
    return res



# res_150kV = max_power_transfer(150e3)
# res_300kV = max_power_transfer(300e3)
# res_400kV = max_power_transfer(400e3)


# plt.figure()
# plt.plot(res_150kV["length"], res_150kV["power"], label="150 kV")
# # plt.plot(res_300kV["length"], res_300kV["power"], label="300 kV")
# plt.plot(res_400kV["length"], res_400kV["power"], label="400 kV")
# plt.xlabel("Cable length [km]")
# plt.ylabel("Transmission active power [MW]")
# plt.title("Transmission capacity vs cable length")
# plt.legend()
# plt.grid()
# plt.show()


######### Test plotting ##########
res_150kV = max_power_transfer_2(150e3)
res_220kV = max_power_transfer_2(220e3)
res_400kV = max_power_transfer_2(400e3)

res_150kV_comp = max_power_transfer(150e3)
res_220kV_comp = max_power_transfer(220e3)
res_400kV_comp = max_power_transfer(400e3)

plt.figure()
plt.plot(res_150kV["length"], res_150kV["power"], label="150 kV", color='tab:blue')
# plt.plot(res_220kV["length"], res_220kV["power"], label="220 kV", color='green')
plt.plot(res_400kV["length"], res_400kV["power"], label="400 kV", color='tab:green')

plt.plot(res_150kV_comp["length"], res_150kV_comp["power"], label="150 kV comp", color='tab:blue', linestyle='--')
# plt.plot(res_220kV["length"], res_220kV["power"], label="220 kV comp", color='green', linestyle='--')
plt.plot(res_400kV_comp["length"], res_400kV_comp["power"], label="400 kV comp", color='tab:green', linestyle='--')

plt.xlabel("Cable length [km]")
plt.ylabel("Transmission active power [MW]")
plt.title("Transmission capacity vs cable length")
plt.legend()
plt.grid()
# plt.show()
plt.savefig("Task2.png")




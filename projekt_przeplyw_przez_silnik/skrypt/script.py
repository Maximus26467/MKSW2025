import cantera as ct
import numpy as np

# Stałe silnika
gamma = 1.4
R = 287.0


# Funkcja do obliczeń parametrów całkowitych z warunków statycznych i liczby Macha
def total_conditions(T_static, p_static, M, gamma=gamma):
    Tt = T_static * (1 + (gamma - 1) / 2 * M**2)
    pt = p_static * (Tt / T_static)**(gamma / (gamma - 1))
    return Tt, pt

# Model silnika turboodrzutowego
class Turbojet:
    def __init__(self, inlet_pressure_ratio, compressor_pressure_ratio,
                 fuel_air_ratio, turbine_pressure_ratio):
        self.pi_d = inlet_pressure_ratio       # Stopień sprężania w dolocie
        self.pi_c = compressor_pressure_ratio  # Stopień sprężania sprężarki
        self.f = fuel_air_ratio                # Stosunek masy paliwa do masy powietrza
        self.pi_t = turbine_pressure_ratio    # Stopień rozprężania turbiny
        self.gas = ct.Solution('gri30.yaml')    # W canterze będziemy korzystać z GRI-Mech 3.0 (gri30.yaml)

    # po wlocie idealnym
    def inlet(self, T0, p0, M0):
        Tt0, pt0 = total_conditions(T0, p0, M0)
        Tt2 = Tt0
        pt2 = pt0 * self.pi_d
        return Tt2, pt2

    # za sprężarką idealną
    def compressor(self, Tt2, pt2):
        pt3 = pt2 * self.pi_c
        Tt3 = Tt2 * (self.pi_c)**((gamma - 1) / gamma)
        return Tt3, pt3

# proces spalania w komorze (paliwem jest metan)
    def combustor(self, Tt3, pt3):
        gas = ct.Solution('gri30.yaml')
        gas.TP = Tt3, pt3

        # Liczba moli powietrza (1 kg) w kmol
        mw_air = gas.mean_molecular_weight
        n_air = 1.0 / mw_air

        # Skład powietrza: 21% O2, 79% N2 molowo
        n_O2 = 0.21 * n_air
        n_N2 = 0.79 * n_air

        # Liczba moli metanu na podstawie masy paliwa f [kg]
        idx_CH4 = gas.species_index('CH4')
        mw_CH4 = gas.molecular_weights[idx_CH4]
        n_CH4 = self.f / mw_CH4

        # Tworzenie mieszaniny molowej
        mixture = {'O2': n_O2, 'N2': n_N2, 'CH4': n_CH4}
        total = sum(mixture.values())
        X = {sp: mixture[sp] / total for sp in mixture}
        gas.X = X

        # Adiabatyczne spalanie pod stałym ciśnieniem
        gas.equilibrate('HP')
        Tt4 = gas.T
        pt4 = pt3
        return Tt4, pt4, gas

    # za turbiną idealną
    def turbine(self, Tt4, pt4):
        pt5 = pt4 / self.pi_t
        Tt5 = Tt4 * (1 / self.pi_t)**((gamma - 1) / gamma)
        return Tt5, pt5

     # za wylotem idealnej dyszy zbieżno rozbieżnej
    def outlet(self, Tt5, pt5, p0):
        Tt0, pt0=total_conditions(T0, p0, M0)
        pi_cr=(2/(gamma+1))**(gamma/(gamma-1))
        pi_n=pt0/pt5
        pt9=pt0
        V9=0.95*((2*gamma*R*Tt5*(1-pi_n**((gamma-1)/gamma)))/(gamma-1))**(1/2)
        Tt9=Tt5-((V9/0.95)**2)*(gamma-1)/(2*gamma*R)
        return Tt9, pt9, V9

# uruchamiamy program
if __name__ == '__main__':
    # dane wejściowe
    T0 = 230.0          # K
    p0 = 30700.0        # Pa
    M0 = 2           # Mach
    pi_d = 0.95         # Stopień sprężania na wlocie
    pi_c = 15.0         # Stopień sprężania sprężarki
    f = 0.2            # Względny wydatek paliwa (kg paliwa / kg powietrza)
    pi_t = 7.3          # Stopień rozprężania turbiny

    engine = Turbojet(pi_d, pi_c, f, pi_t)

    # 1. Za wlotem
    Tt2, pt2 = engine.inlet(T0, p0, M0)
    print(f"Za wlotem: Tt2 = {Tt2:.2f} K, pt2 = {pt2:.2f} Pa")

    # 2. Za sprężarką
    Tt3, pt3 = engine.compressor(Tt2, pt2)
    print(f"Po sprężarce: Tt3 = {Tt3:.2f} K, pt3 = {pt3:.2f} Pa")

    # 3. W komorze spalania
    Tt4, pt4, gas_out = engine.combustor(Tt3, pt3)
    print(f"Po spalaniu: Tt4 = {Tt4:.2f} K, pt4 = {pt4:.2f} Pa")

    # 4. Za turbiną
    Tt5, pt5 = engine.turbine(Tt4, pt4)
    print(f"Po turbinie: Tt5 = {Tt5:.2f} K, pt5 = {pt5:.2f} Pa")

    # 4. Za wylotem
    Tt9, pt9, V9 = engine.outlet(Tt5, pt5, p0)
    print(f"Za wylotem: Tt9 = {Tt9:.2f} K, pt9 = {pt9:.2f} Pa, V9 = {V9:.2f} m/s")

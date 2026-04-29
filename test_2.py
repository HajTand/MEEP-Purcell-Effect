import math
import matplotlib.pyplot as plt
import numpy as np

import meep as mp

"количество микселей на микрометров"
"размер ячейки 5 нм = 1/200"
resolution = 200  # pixels/μm
"показатель преломления среды"
n = 2.4
"диэлектрическая проницаемость среды"
epsilon = n**2

"длина волны и частота диполя в вакууме"
dipole_vac_wavelen = 1.0    #μm
dipole_vac_freq = 1/dipole_vac_wavelen #MHz

"dpml - граничное условие, чем больше, тем лучше поглощает (не отражает)"
"d - значит слой в микрометрах"
dpml = 1.0

"Достаточно широкий, чтобы импульс быстро затух (симуляция быстрая)"
"Достаточно узкий, чтобы не возбуждать далекие частоты (точность)"
"Ширина испульса диполя"
df = 0.2 * dipole_vac_freq
nfreq = 1


# расчет без зеркал, просто однородная среда

"L - размер среды"
L = 10.0
"sz - размер рабочей области целого вместе с pml"
sz = L + 2.0 * dpml
cell_size = mp.Vector3(z=sz)
"размер ячейки"
pml_layers = [mp.PML(dpml)]
"количество слоев pml"

"задаем источник и его параметры, центр геометрический"
"ширины df источник"
sources = [ 
    mp.Source(
        mp.GaussianSource(dipole_vac_freq, fwidth=df),
        component=mp.Ex,
        center=mp.Vector3(z=0),
    )
]

"задаем симуляцию, dim - размерность моделирования"
"размер ячейки - cell_size, исчтоник - sources, разрешение - resolution"
"boundary layers - граничные слои"
sim = mp.Simulation(
    cell_size=cell_size,
    boundary_layers=pml_layers,
    sources=sources,
    dimensions=1,
    resolution=resolution,
    default_material=mp.Medium(index=n)
)

"запуск симуляции, задаем когда закончится симуляция"
"каждые 50 фемтосекунд после того как закончился импульс света"
"компонента Ex должно стать 1e-9"
sim.run(mp.dft_ldos(dipole_vac_freq, df, nfreq),
        until_after_sources=mp.stop_when_fields_decayed(50, mp.Ex, mp.Vector3(), 1e-9))

ldos_bulk = sim.ldos_data[0]

sim.reset_meep()

# рассчет с резонатором
f_P = []
Norm_Thickness = np.arange(0.2,3.2,0.2)

for value in Norm_Thickness:
    sz = ((value * dipole_vac_wavelen) / n)
    cell_size = mp.Vector3(z = sz)
    sources = [ 
        mp.Source(
            mp.GaussianSource(dipole_vac_freq, fwidth=df),
            component=mp.Ex,
            center=mp.Vector3(z=0),
        )
    ]
    
    sim = mp.Simulation(
       cell_size=cell_size,
       sources=sources,
       dimensions=1,
       resolution=resolution,
       default_material=mp.Medium(index=n)
       )
       # boundary_layers НЕТ! потому что границы = PEC зеркала
    
    sim.run(mp.dft_ldos(dipole_vac_freq, df, nfreq),
            until_after_sources=1000)

    ldos_res = sim.ldos_data[0]
    f_P.append(ldos_res/ldos_bulk)
    sim.reset_meep()

#теория
def purcell_theory(x):
    """
    Аналитическая формула из IEEE J. Quantum Electronics, Vol. 34, pp. 71-76 (1998)
    x = nL/λ  — нормированная толщина резонатора
    """
    m = np.floor(x + 0.5)  # номер моды
    term1 = (3 * m) / (4 * x)
    term2 = (4 * m**3 - m) / (16 * x**3)
    return term1 + term2

# Рассчитываем теоретические значения для тех же точек
purcell_theory_vals = [purcell_theory(x) for x in Norm_Thickness] 
   
print()
plt.plot(Norm_Thickness,f_P)
plt.plot(Norm_Thickness,purcell_theory_vals)

print(ldos_bulk)



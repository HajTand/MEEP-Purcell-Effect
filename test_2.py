import math
import matplotlib.pyplot as plt
import numpy as np

import meep as mp

"количество микселей на микрометров"
"размер ячейки 5 нм = 1/200"
resolution = 200  # pixels/μm
"показатель преломления среды"
n = 2.4

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

"sz - размер рабочей области целого вместе с pml"
sz = 10.0 + 2.0 * dpml
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
        center=mp.Vector3(z=-0.5 * sz + dpml),
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
)

sim.run(mp.dft_ldos(dipole_vac_freq, df, 1),
            until_after_sources=mp.stop_when_fields_decayed(50, mp.Ex, mp.Vector3(), 1e-6))
ldos_sim0 = sim.ldos_data[0]
print(ldos_sim0)
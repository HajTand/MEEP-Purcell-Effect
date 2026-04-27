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
sz = 10 + 2 * dpml
cell_size = mp.Vector3(z=sz)


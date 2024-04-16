import lifesim

# ----- README
# - This file is used to execute an experiment.
# - here we attempt to find the best configuration of the emma x array.
# - we start with the emma x array, as there have been experiments done on this before,
# and so we have the opportunity to validate the results
# - the only changes are made using the Options file
# First we apply changes through a grid search, then through an algorithm
# NOPE sometimes you gotta run before you can walk

# imports
import time
import pygmo as pg
import numpy as np
import lifesim.util.figures as fg
import os
import itertools
from tqdm import tqdm
import matplotlib.pyplot as plt
from lifesim.util.normalizer import baseline_distances
# --------------- modules
bus = lifesim.Bus()
instrument = lifesim.Instrument(name='inst')
bus.add_module(instrument)
transm = lifesim.TransmissionMap(name='transm')
bus.add_module(transm)
exo = lifesim.PhotonNoiseExozodi(name='exo')
bus.add_module(exo)
local = lifesim.PhotonNoiseLocalzodi(name='local')
bus.add_module(local)
star = lifesim.PhotonNoiseStar(name='star')
bus.add_module(star)
ins_shot = lifesim.PhotonNoiseThermalInst(name='ins_shot')
bus.add_module(ins_shot)
# --------------- connect all modules
bus.connect(('inst', 'transm'))
bus.connect(('inst', 'exo'))
bus.connect(('inst', 'local'))
bus.connect(('inst', 'star'))
bus.connect(('inst', 'ins_shot'))
bus.connect(('star', 'transm'))
# --------------- options
bus.data.options.set_scenario()
# --------------- import test planet
Earth_clone = np.array([1.0, 365, 1.0, 0.0, 0, 0, 0, 0., 0.29,
                        0, 0, 3.0, 1.0, 1.0, 0.1, 0, 0, 0,
                        285, 1.0, 0, 5778, 2.5, 0.0, 0., 0, 0, 0, 0]).reshape((1, 29))

bus.data.catalog_add_planet(Earth_clone)
bus.data.catalog.loc[0, 'lon'] = 135.           # this is the proper format for indexing dataframes
bus.data.catalog.loc[0, 'lat'] = 45.
bus.data.catalog.loc[0, 'angsep'] = bus.data.catalog.loc[0, 'sep_p']/bus.data.catalog.loc[0, 'distance_s']
print(bus.data.catalog.iloc[0])
# ----------------- ouput file
root = 'C:/Users/rogie/Desktop/Afstudeerproject/09. Experiments/02. Output/single/'
run_name = 's2'
save_bool = True
# ---------------  INPUTS
grid = [0.1, 0.25, 0.5, 0.75, 1.0]
configs = ['TTN', 'emma x', 'angel woolf', 'diamond', 'kite', 'pentagon']
preset_configs = [[2, 3, 'TTN'],
                   [1, 4, 'angel woolf'],
                   [1, 4, 'emma x'],
                   [1, 4, 'diamond'],
                   [2, 4, 'kite'],
                   [2, 5, 'pentagon']]
baselines = np.array(grid)*2.
ratios = np.array(grid)*6.
areas = np.array(grid)*4.
temperatures = np.array(grid)*80.
fom_array = np.zeros((6, 5, 5, 5, 5))
# loop outputs
for do, preset in enumerate(tqdm(preset_configs)):
    print('Ladies and gentlemen, config nr ', do)
    for re, baseline in enumerate(baselines):
        for mi, ratio in enumerate(ratios):
            for la, area in enumerate(areas):
                for ti, temperature in enumerate(temperatures):
                    bus.data.options.set_manual(technique=preset[0],
                                                n_apertures = preset[1],
                                                array_config = preset[2],
                                                baseline=baseline,
                                                ratio=ratio,
                                                area=area,
                                                instrument_temperature=temperature)
                    instrument.apply_options()
                    mirror_radius = (area / preset[1] / np.pi) ** 0.5
                    dimension_check, _ = baseline_distances(bus.data.inst['configuration'])
                    if mirror_radius >= dimension_check * 0.5:
                        print('---- rejected instrument')
                        continue
                    instrument.get_snr()
                    fom = bus.data.catalog.snr_1h
                    fom_array[do, re, mi, la, ti] = fom
print('done')

if save_bool:
    import pickle
    with open(root + run_name+ '_earth_snrs.pkl', 'wb') as f:
        pickle.dump(fom_array, f)


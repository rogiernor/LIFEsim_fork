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
# --------------- import test pop
root = 'C:/Users/rogie/Desktop/Afstudeerproject/03. Code/01. LIFEsim-master/'
planet_path = os.path.join(root, 'docs')
bus.data.import_catalog(
    input_path=planet_path+'/FOM_ppop.hdf5')
print('done',  len(bus.data.catalog))
bus.data.catalog_remove_distance(stype=0, mode='larger', dist=0.)  # remove all A stars
bus.data.catalog_remove_distance(stype=4, mode='larger', dist=5.)  # remove M stars > 5pc
bus.data.catalog_remove_distance(stype=3, mode='larger', dist=5.)  # remove K stars > 5pc
bus.data.catalog_remove_distance(stype=2, mode='larger', dist=5.)  # remove G stars > 5pc
bus.data.catalog_remove_distance(stype=1, mode='larger', dist=5.)  # remove F stars > 5pc
# --------------- Select population of 558 planets
# bus.data.catalog.info()
ppop_pop_size = 500
ppop_parser = len(bus.data.catalog)//ppop_pop_size
# print('the amount of skips', ppop_parser)
bus.data.catalog = bus.data.catalog.loc[::ppop_parser]
# ----------------- ouput file
root2 = 'C:/Users/rogie/Desktop/Afstudeerproject/09. Experiments/02. Output/Emma circle/'
run_name = 's5'
save_bool = True
preset_configs = [1, 4, 'diamond']
area = 0.75
radius = 4
temperature = 50
# --- iunitialize
baselines = np.arange(1.3, 6.6, 0.1)
# baselines = np.arange(0.25, 10.0, 0.25)
# loop outputs
for i, baseline in enumerate(baselines):
    ratio = np.sqrt(radius**2/(0.25*baseline**2)-1)
    bus.data.options.set_manual(technique=preset_configs[0],
                                n_apertures = preset_configs[1],
                                array_config = preset_configs[2],
                                baseline=baseline,
                                ratio=ratio,
                                area=area,
                                instrument_temperature=temperature)
    instrument.apply_options()
    mirror_radius = (area / preset_configs[1] / np.pi) ** 0.5
    dimension_check, _ = baseline_distances(bus.data.inst['configuration'])
    if mirror_radius >= dimension_check * 0.5:
        print('---- rejected instrument')
        continue
    instrument.get_snr()
    fom = bus.data.catalog.snr_1h
    if save_bool:
        test_name = 'diamond_33bl_075A_'+str(baseline*100)
        output_path = os.path.join(root2, test_name)
        bus.data.export_catalog_txt(output_path=output_path + '.txt')
        bus.data.export_catalog(output_path=output_path + '.hdf5')
        bus.data.export_instrument(output_path=output_path)
print('done')


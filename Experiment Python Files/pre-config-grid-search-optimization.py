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
# --------------- import test pop #TODO export this FOM to a method in DATA
root = 'C:/Users/rogie/Desktop/Afstudeerproject/09. Experiments/02. Output/'

planet_path = "C:/Users/rogie/Desktop/Afstudeerproject/03. Code/01. LIFEsim-master/docs/FOM_ppop.hdf5"
bus.data.import_catalog(
    input_path=planet_path)
print('done',  len(bus.data.catalog))
bus.data.catalog_remove_distance(stype=0, mode='larger', dist=0.)  # remove all A stars
bus.data.catalog_remove_distance(stype=4, mode='larger', dist=5.)  # remove M stars > 10pc
bus.data.catalog_remove_distance(stype=3, mode='larger', dist=5.)  # remove K stars > 10pc
bus.data.catalog_remove_distance(stype=2, mode='larger', dist=5.)  # remove G stars > 10pc
bus.data.catalog_remove_distance(stype=1, mode='larger', dist=5.)  # remove F stars > 10pc
# --------------- Select population of 50 planets
bus.data.catalog.info()
ppop_pop_size = 500
ppop_parser = len(bus.data.catalog)//ppop_pop_size
# print('the amount of skips', ppop_parser)
bus.data.catalog = bus.data.catalog.loc[::ppop_parser]
# -------------TODO All simulation settings
# experiment_settings = [[2, 3, 'TTN'],
#                        [1, 4, 'angel woolf'],
#                        [1, 4, 'emma x'],
#                        [1, 4, 'diamond'],
#                        [2, 4, 'kite'],
#                        [2, 5, 'pentagon']]
experiment_settings = [[1, 4, 'angel woolf']]
#
parameters = [2., 6., 4., 80]
combinations = list(itertools.product([0.1, 0.25, 0.5, 0.75, 1.0], repeat=4))
save_bool = True
#--------------
print('start at', time.perf_counter())
for n, setting in enumerate(experiment_settings):
    scenario = setting[2]
    bus.data.options.set_manual(technique=setting[0])
    bus.data.options.set_manual(n_apertures=setting[1])
    bus.data.options.set_manual(array_config=setting[2])

    if save_bool:
        path = root+'design space output/grid_search_30_1_fom_500_constrained'
        directory = os.path.join(path, scenario+'2')
        scenario_file = os.mkdir(directory)
    print('this run simulates an instrument which uses technique '+str(setting[0])+' with '+str(setting[1])+' apertures'
          ' in a(n) '+setting[2]+' config to loop over 500 close-by planets. Now lets create this instrument')
    # fg.configuration_render(bus.data.inst['configuration'], bus.data.inst['diameters'], save=False,
    #                         save_name=scenario+' configuration a_e_e')
    # print('Total Array: ', bus.data.options.array)
    for m, combination in enumerate(combinations):
        baseline = parameters[0]*combination[0]
        ratio = parameters[1]*combination[1]
        area = parameters[2]*combination[2]
        temperature = parameters[3]*combination[3]
        print('which yields for '+str(m)+' boundary values baseline: ', baseline, ', ratio: ', ratio,
              ', area: ', area, 'temperature', temperature)
        bus.data.options.set_manual(baseline=baseline,
                                    ratio=ratio,
                                    area=area,
                                    instrument_temperature=temperature)
        # create instrument
        instrument.apply_options()
        # --- check instrument
        mirror_radius = (area / bus.data.options.array['n_apertures'] / np.pi) ** 0.5
        dimension_check, _ = baseline_distances(bus.data.inst['configuration'])
        if mirror_radius >= dimension_check*0.5:
            print('---- rejected instrument')
            continue
        # fg.configuration_render(bus.data.inst['configuration'], bus.data.inst['diameters'], save=False,
        #                         save_name=scenario + ' configuration a_e_e')
        # ---- verbose instrument
        # print('area of instrument ', bus.data.inst['telescope_area'])
        # print('instrument after generation', bus.data.inst['configuration'])
        # ------------ Run calculations
        # instrument.get_snr(verbose=False)
        # # # ---------- Saving the Results ----------
        # print('which results in snrs ', print(bus.data.catalog['snr_1h']))
        # if save_bool:
        #     directory = 'output/Grand comparison/equal area (2 m total)'
        #     test_name = scenario+'_earth_clone_distances_1_10_pc'
        #     output_path = os.path.join(root, directory, test_name)
        #     bus.data.export_catalog_txt(output_path=output_path+'.txt')
        #     bus.data.export_catalog(output_path=output_path+'.hdf5')
        #     bus.data.export_instrument(output_path=output_path)
        # -------- Figure configuration
        # fg.configuration_render(bus.data.inst['configuration'], bus.data.inst['diameters'], save=False,
        #                         save_name=scenario+' configuration a_e_e')
        # wl_slab = 0
        # instrument.run_socket(s_name='transmission', method='plot_psf', slab=wl_slab, save=False,
        #                       save_name=scenario+' psf plot')
        # -------- Figure of all
        instrument.get_snr(verbose=False)
        print('~done!~')
        if save_bool:
            test_name = '557_5_pc_planets_'+str(baseline)+'_'+str(ratio)+'_'+str(area)+'_'+str(temperature)
            output_path = os.path.join(directory, test_name)
            bus.data.export_catalog_txt(output_path=output_path+'.txt')
            bus.data.export_catalog(output_path=output_path+'.hdf5')
            bus.data.export_instrument(output_path=output_path)
    print('~~~next config~~~', time.perf_counter())
print('end of render @', time.perf_counter())

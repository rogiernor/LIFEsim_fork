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
bus.connect(('inst', 'exo'))  # TODO How the fuck does this happen
bus.connect(('inst', 'local'))
bus.connect(('inst', 'star'))
bus.connect(('inst', 'ins_shot'))
bus.connect(('star', 'transm'))
# --------------- options
bus.data.options.set_scenario()
# --------------- import test pop #TODO export this FOM to a method in DATA
root = 'C:/Users/rogie/Desktop/Afstudeerproject/03. Code/01. LIFEsim-master/'
planet_path = os.path.join(root, 'docs')
bus.data.import_catalog(
    input_path=planet_path+'/FOM_ppop.hdf5')
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
print('the amount of skips', ppop_parser)
bus.data.catalog = bus.data.catalog.loc[::ppop_parser]
print('catlog, dude', bus.data.catalog)

# alternative fom
# bus.data.catalog = bus.data.catalog.loc[bus.data.catalog['nuniverse'] == 0]
# bus.data.catalog = bus.data.catalog.loc[bus.data.catalog['habitable'] == True]
print(' FOM', bus.data.catalog)
# -------------TODO All simulation settings
# experiment_settings = [[2, 3, 'TTN'],
#                        [1, 4, 'angel woolf'],
#                        [1, 4, 'emma x'],
#                        [1, 4, 'diamond'],
#                        [2, 4, 'kite'],
#                        [2, 5, 'pentagon']]

experiment_settings = [[2, 3, 'TTN'],[1, 4, 'emma x']]


# -------------------
class X_optimization:
    def __init__(self, dim):
        self.dim = dim

    def fitness(self, x):
        area = x[2]
        bus.data.options.set_manual(n_apertures=experiment_setting[1],
                                    technique=experiment_setting[0],
                                    array_config=experiment_setting[2],
                                    baseline=x[0],
                                    ratio=x[1],
                                    area=area,
                                    instrument_temperature=50)

        instrument.apply_options()
        # fg.configuration_render(bus.data.inst['configuration'],bus.data.inst['diameters'])
        instrument.get_snr(verbose=False)
        # Aperture constraint
        mirror_radius = (area / experiment_setting[1] / np.pi) ** 0.5
        dimension_check, _ = baseline_distances(bus.data.inst['configuration'])
        # radius constraint
        radii_of_config = np.linalg.norm(bus.data.inst['configuration'], axis=-1)
        fom2 = max(radii_of_config)
        print('radii', radii_of_config)
        constraint = 0
        if mirror_radius >= dimension_check*0.5:
            print('---- rejected instrument: overlap')
            constraint = 1000
        if fom2 >= 5.:
            print('---- rejected instrument: too large')
            constraint = 1000
        # if fom2 <= 3.1:
        #     print('---- rejected instrument: too small')
        #     constraint = 1000

        fom = len(bus.data.catalog[bus.data.catalog.snr_1h < 5.])
        # fom = bus.data.catalog[bus.data.catalog.snr_1h < 5.]
        # print(fom[['radius_p', 'sep_p', 'distance_s', 'temp_p', 'habitable', 'snr_1h']]) #, 'dista
        # print(len(fom))
        # fom = len(bus.data.catalog[(bus.data.catalog.radius_p < 2) & (bus.data.catalog.snr_1h < 5)])
        # print('fom2', fom2)
        # fom = len(bus.data.catalog[(bus.data.catalog.radius_p < 2) & (bus.data.catalog.snr_1h < 5)])
        # print('fom', fom)
        # print('fom2', fom2)
        return [fom, x[2], constraint]

    def get_bounds(self):
        return [0.1, 0.1, 0.1], [10.0, 2.0, 4.0]

    def get_nobj(self):
        return 3

    def get_name(self):
        return "X optimization"

    def get_extra_info(self):
        return "\tDimensions: " + str(self.dim)

bus.data.options.set_manual(random_seed=1)
save_zone = 'C:/Users/rogie/Desktop/Afstudeerproject/05. Figures (made by me)/02. Thesis/03. Result figures/03. MO optimization fronts/'
root2 = 'C:/Users/rogie/Desktop/Afstudeerproject/09. Experiments/02. Output/MO pops/'
mo_pops = ['mv57','mv59']
ylabel = 'Max Area [m]'
xlabel = 'Total non-detections [-]'
# configs = [0, 1]
save_bool = True
for i, run_number in enumerate(mo_pops):
    print(run_number)
    experiment_setting = experiment_settings[i]
    udp = pg.problem(X_optimization(3))
    algo = pg.algorithm(pg.maco(gen=1, ker=9, seed=bus.data.options.other['random_seed']))
    pop = pg.population(udp, 20, seed=bus.data.options.other['random_seed'])
    pop1 = algo.evolve(pop)
    pop2 = algo.evolve(pop1)
    pop3 = algo.evolve(pop2)
    pop4 = algo.evolve(pop3)
    pop5 = algo.evolve(pop4)
    print(pop5)
    im2 = pg.plot_non_dominated_fronts(pop5.get_f())
    im2.set(xlabel=xlabel, ylabel=ylabel)
    plt.savefig(save_zone + run_number + 'gen_3')
    plt.show()
    if save_bool:
        import pickle
        with open(root2 + run_number + '_eval_pop.pkl', 'wb') as f:
            pickle.dump(pop5, f)
print('DONE')
# # ----------- gen 1
# pop1 = algo.evolve(pop)
# print(pop1)
# im1 = pg.plot_non_dominated_fronts(pop1.get_f())
# im1.set(xlabel=xlabel, ylabel=ylabel)
# plt.savefig(save_zone+run_number+'gen_1')
# plt.show()
# # ---------- gen 2
# pop2 = algo.evolve(pop1)
# im2 = pg.plot_non_dominated_fronts(pop2.get_f())
# im2.set(xlabel=xlabel, ylabel=ylabel)
# plt.savefig(save_zone+run_number+'gen_2')
# plt.show()
# print(pop2)
# # ---------- gen 3
# pop3 = algo.evolve(pop2)
# im2 = pg.plot_non_dominated_fronts(pop3.get_f())
# im2.set(xlabel=xlabel, ylabel=ylabel)
# plt.savefig(save_zone+run_number+'gen_3')
# plt.show()
# print(pop3)
# # ---------- gen 3
#
#
# ndf, dl, dc, ndr = pg.fast_non_dominated_sorting(pop3.get_f())
# print('front 1', ndf)
# print('domination list', dl)
# print('domination counter', dc)
# print('non dominated ranks', ndr)

#
# if save_bool:
#         import pickle
#         with open(root2 + run_number + '_eval_pop.pkl', 'wb') as f:
#             pickle.dump(pop3, f)
#
#
#         with open(root2 + run_number + '_eval_pop.pkl', 'rb') as f:
#             pop_loaded = pickle.load(f)
#
# print('retrieve',pop_loaded.get_x())


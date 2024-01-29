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
bus.data.catalog_remove_distance(stype=4, mode='larger', dist=5.)  # remove M stars > 10pc
bus.data.catalog_remove_distance(stype=3, mode='larger', dist=5.)  # remove K stars > 10pc
bus.data.catalog_remove_distance(stype=2, mode='larger', dist=5.)  # remove G stars > 10pc
bus.data.catalog_remove_distance(stype=1, mode='larger', dist=5.)  # remove F stars > 10pc
# --------------- Select population of 50 planets
# bus.data.catalog.info()
ppop_pop_size = 500
ppop_parser = len(bus.data.catalog)//ppop_pop_size
# print('the amount of skips', ppop_parser)
bus.data.catalog = bus.data.catalog.loc[::ppop_parser]
print(len(bus.data.catalog))
# ---------------- Run SNR
print('from modules', time.perf_counter())
# bus.data.options.set_manual(baseline=2.0, n_apertures=4, technique=1, array_config='emma x')
# instrument.apply_options()
# instrument.get_snr(verbose=False)
# experiment_settings = [[2, 3, 'TTN'],
#                        [1, 4, 'angel woolf'],
#                        [1, 4, 'emma x'],
#                        [1, 4, 'diamond'],
#                        [2, 4, 'kite'],
#                        [2, 5, 'pentagon']]

class X_optimization:
    def __init__(self, dim):
        self.dim = dim

    def fitness(self, x):
        print('genome', x)
        area = 2
        experiment_settings = [[1, 4, 'emma x'],
                               [1, 4, 'diamond'],
                               [2, 4, 'kite']]

        max_rad = 6
        # ratio = np.sqrt(max_rad**2/(0.25*x[0]**2)-1)
        ratio = x[1]
        print('Config check:', experiment_settings[config_index][2])
        bus.data.options.set_manual(n_apertures=experiment_settings[config_index][1],
                                    technique=experiment_settings[config_index][0],
                                    array_config=experiment_settings[config_index][2],
                                    baseline=x[0],
                                    ratio=ratio,
                                    area=area,
                                    instrument_temperature=50)
        instrument.apply_options()
        # ---- CONSTRAINTS
        # --- check instrument
        mirror_radius = (area / self.dim / np.pi) ** 0.5 # n_ap IS CALCULATED AT int(3)! So radii nare always slightly larger
        dimension_check, _ = baseline_distances(bus.data.inst['configuration'])
        constraint = 0
        if mirror_radius >= dimension_check*0.5:
            print('---- rejected instrument')
            constraint = 1
        # ---- check max baseline, less than 10
        # in_constraint = x[0] * x[1] - 10
        radii_of_config = np.linalg.norm(bus.data.inst['configuration'], axis=-1)
        constraint2 = max(radii_of_config)-max_rad
        print('constraint values, number 1:', constraint, 'constraint value 2:', radii_of_config)
        # ---- FITNESS
        instrument.get_snr(verbose=False)
        # fom = len(bus.data.catalog[bus.data.catalog.snr_1h > 5])
        # reduce to only < 2 radius planetsbus.data.catalog.snr_1h < 5
        # fom = len(bus.data.catalog[(bus.data.catalog.radius_p < 2) & (bus.data.catalog.snr_1h < 5)])
        # reduce to only habitable
        # fom = len(bus.data.catalog[(bus.data.catalog.habitable == True) & (bus.data.catalog.snr_1h < 5)])
        fom = len(bus.data.catalog[bus.data.catalog.snr_1h < 5])
        # print(fom)
        return [fom, constraint, constraint2]

    def get_nec(self):
        return 1
    #
    def get_nic(self):
        return 1

    # def get_nix(self):
    #     return 1

    def get_bounds(self):
        return [0.1, 0.1], [10.0, 2.0]

    def get_name(self):
        return "X optimization"

    def get_extra_info(self):
        return "\tDimensions: " + str(self.dim)

print('from to end of run', time.perf_counter())
root2 = 'C:/Users/rogie/Desktop/Afstudeerproject/09. Experiments/02. Output/winner pops/'
bus.data.options.set_manual(random_seed=1)
# initialize optimization
run_numbers = ['74', '75', '76']
indexes = [0,1,2]
save_bool = True
for i in indexes:
    run_number = run_numbers[i]
    config_index = i
    print('We are doing config:', i)
    prob = pg.problem(X_optimization(2))
    print(prob)
    algo = pg.algorithm(pg.gaco(gen=1, ker=9, seed=bus.data.options.other['random_seed']))
    pop = pg.population(prob, 20, seed=bus.data.options.other['random_seed'])
    pop1 = algo.evolve(pop)
    # print(pop)
    pop2 = algo.evolve(pop1)
    # print(pop)
    pop3 = algo.evolve(pop2)
    # print(pop)
    print('Winner pop for config', i, 'is', pop3)
    if save_bool:
        import pickle
        with open(root2 + run_number + '_eval_pop.pkl', 'wb') as f:
            pickle.dump(pop3, f)
print('DONE')





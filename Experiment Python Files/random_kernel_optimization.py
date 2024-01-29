import lifesim
# ----- README
# - This file is used to execute an experiment.
#

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
bus.data.catalog_remove_distance(stype=4, mode='larger', dist=5.)  # remove M stars > 5pc
bus.data.catalog_remove_distance(stype=3, mode='larger', dist=5.)  # remove K stars > 5pc
bus.data.catalog_remove_distance(stype=2, mode='larger', dist=5.)  # remove G stars > 5pc
bus.data.catalog_remove_distance(stype=1, mode='larger', dist=5.)  # remove F stars > 5pc
# --------------- Select population of 50 planets
# bus.data.catalog.info()
ppop_pop_size = 500
ppop_parser = len(bus.data.catalog)//ppop_pop_size
# print('the amount of skips', ppop_parser)
bus.data.catalog = bus.data.catalog.loc[::ppop_parser]
# ---------------- Random Kernel Optimization
print('from modules', time.perf_counter())


class random_optimization:
    def __init__(self, dim, loose: int = 0):
        print(dim)
        self.dim = dim
        self.loose = loose

    def fitness(self, x):
        r = [4.0]*(self.dim)
        if self.loose >= 1:
            r[-self.loose:] = x[-self.loose:]
        area = 0.75
        geometry_list = []
        for i in list(range(self.dim)):
            aperture_position = [r[i] * np.cos(x[i]), r[i] * np.sin(x[i])]
            geometry_list.append(aperture_position)
        geometry = np.array(geometry_list)
        bus.data.options.set_manual(n_apertures=self.dim, area=area, technique=2, array_config=geometry)
        instrument.apply_options()
        # limit = 6
        # fg.configuration_render(bus.data.inst['configuration'], bus.data.inst['diameters'], dimensions=limit)
        # --- check instrument
        mirror_radius = (area / self.dim / np.pi) ** 0.5
        dimension_check, _ = baseline_distances(bus.data.inst['configuration'])
        constraint = 0
        if mirror_radius >= dimension_check*0.5:
            print('---- rejected instrument')
            constraint = 1
        instrument.get_snr(verbose=False)
        fom = len(bus.data.catalog[bus.data.catalog.snr_1h < 5])
        return [fom, constraint]

    def get_bounds(self):
        bottom = [0.0] * self.dim
        top = [2 * np.pi] * self.dim
        if self.loose >= 1:
            bottom.extend([0.0] * self.loose)
            top.extend([8.0] * self.loose)
        print(bottom,top)
        return bottom, top

    def get_nec(self):
        return 1

    def get_name(self):
        return "X optimization"

    def get_extra_info(self):
        return "\tDimensions: " + str(self.dim)


print('from to end of run', time.perf_counter())
root2 = 'C:/Users/rogie/Desktop/Afstudeerproject/09. Experiments/02. Output/winner pops/'
save_bool = True
dimmadome = [3,4,5]
run_numbers = ['r35', 'r36', 'r37']
bus.data.options.set_manual(random_seed=1)
for i, run_number in enumerate(run_numbers):
    prob = pg.problem(udp = random_optimization(dimmadome[i],0))
    # prob = pg.unconstrain(prob=prob, method='death penalty')
    print(prob)
    algo = pg.algorithm(pg.gaco(gen=1, ker=9, seed=bus.data.options.other['random_seed']))
    # algo.set_verbosity(2)
    pop = pg.population(prob, 20, seed=bus.data.options.other['random_seed'])
    # print(pop)
    pop1 = algo.evolve(pop)
    # print(pop)
    pop2 = algo.evolve(pop1)
    # print(pop)
    pop3 = algo.evolve(pop2)
    # print(pop)
    pop4 = algo.evolve(pop3)
    # print(pop)
    pop5 = algo.evolve(pop4)
    # print(pop)
    pop6 = algo.evolve(pop5)
    print(pop6)
    if save_bool:
        import pickle
        with open(root2 + run_number + '_eval_pop.pkl', 'wb') as f:
            pickle.dump(pop6, f)
# im = pg.plot_non_dominated_fronts(pop.get_f())
# im.set(xlabel='Planetary non detection', ylabel='', )
# plt.savefig(save_zone+run_number+'gen_2')
# plt.show()

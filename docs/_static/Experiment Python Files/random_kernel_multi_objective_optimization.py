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
# -----------  zoom everything in
# bus.data.catalog = bus.data.catalog.loc[bus.data.catalog['habitable'] == True]
# print('HABITABLE FOM', bus.data.catalog)
# ---------------- Run SNR
print('from modules', time.perf_counter())




class random_mo_optimization:
    def __init__(self, dim):
        self.dim = dim

    def fitness(self, x):
        r = x[-1]
        area = 2
        n_apertures = self.dim-1
        geometry_list = []
        for i in list(range(self.dim-1)):
            aperture_position = [r * np.cos(x[i]), r * np.sin(x[i])]
            geometry_list.append(aperture_position)
        geometry = np.array(geometry_list)
        bus.data.options.set_manual(n_apertures=n_apertures, area=area, technique=2, array_config=geometry)
        instrument.apply_options()
        # --- check instrument
        mirror_radius = (area / n_apertures / np.pi) ** 0.5
        dimension_check, _ = baseline_distances(bus.data.inst['configuration'])
        instrument.get_snr(verbose=False)
        fom = len(bus.data.catalog[bus.data.catalog.snr_1h < 5.])
        # fom = len(bus.data.catalog[(bus.data.catalog.radius_p < 2) & (bus.data.catalog.snr_1h < 5)])
        constraint = 0
        if mirror_radius >= dimension_check*0.5:
            print('---- rejected instrument')
            constraint = 100
        return [fom, x[-1], constraint]

    def get_nobj(self):
        return 3

    def get_bounds(self):
        bottom = [0.0]*self.dim
        top = [2*np.pi]*self.dim
        top[-1] = 10 # Radius
        return bottom, top


    def get_name(self):
        return "Kernel nulling MO"

    def get_extra_info(self):
        return "\tDimensions: " + str(self.dim)


print('from to end of run', time.perf_counter())
bus.data.options.set_manual(random_seed=1)
# generate algorithm.
root2 = 'C:/Users/rogie/Desktop/Afstudeerproject/09. Experiments/02. Output/MO pops/'
save_zone = 'C:/Users/rogie/Desktop/Afstudeerproject/05. Figures (made by me)/02. Thesis/03. Result figures/03. MO optimization fronts/'
run_numbers = ['rando_single_pop_evolve_twice']
dimensions = [4]
save_bool = True
xlabel = 'Planetary non-detections [-]'
ylabel = 'Radius [m]'

for i, run_number in enumerate(run_numbers):
    # --------- begin
    udp = pg.problem(random_mo_optimization(dimensions[i]))
    algo = pg.algorithm(pg.maco(gen=1, ker=9, seed=bus.data.options.other['random_seed']))
    print(udp)
    # gen ------------- 0
    pop = pg.population(udp, 20, seed=bus.data.options.other['random_seed'])
    pop1 = algo.evolve(pop)
    pop2 = algo.evolve(pop1)
    # pop3 = algo.evolve(pop2)
    # pop4 = algo.evolve(pop3)
    if save_bool:
        import pickle
        with open(root2 + run_number + '_eval_pop.pkl', 'wb') as f:
            pickle.dump(pop2, f)
print('DONE')
# ndf, dl, dc, ndl = pg.fast_non_dominated_sorting(pop.get_f())
# print('front 1', ndf, ' dl ', dl, 'dc', dc, 'ndl', ndl)
# im = pg.plot_non_dominated_fronts(pop.get_f())
# im.set(xlabel='Planetary non-detections [-]', ylabel='nulling baseline [m]')
# plt.show()
# # ----------- gen 1
# pop1 = algo.evolve(pop)
# im1 = pg.plot_non_dominated_fronts(pop1.get_f())
# im1.set(xlabel='Planetary non-detections [-]', ylabel='nulling baseline [m]')
# plt.show()
# # ---------- gen 2
# pop2 = algo.evolve(pop1)
# im2 = pg.plot_non_dominated_fronts(pop2.get_f())
# im2.set(xlabel='Planetary non-detections [-]', ylabel='nulling baseline [m]')
# plt.show()
# # ---------- gen 3
# pop3 = algo.evolve(pop2)
# im2 = pg.plot_non_dominated_fronts(pop3.get_f())
# im2.set(xlabel=xlabel, ylabel=ylabel)
# plt.savefig(save_zone+run_number+'gen_3')
# plt.show()
# print(pop3)



# if save_bool:
#         import pickle
#         with open(root2 + run_number + '_eval_pop.pkl', 'wb') as f:
#             pickle.dump(pop3, f)
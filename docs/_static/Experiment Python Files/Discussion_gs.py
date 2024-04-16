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
experiment_settings = [[1, 4, 'emma x']]
#
what_are_we_varying = '_time_5mr_A0_75m'
save_bool = True
figure_bool = True
#
# x_axis = np.array([0.01, 0.05, 0.1, 0.15, 0.20, 0.25, 0.3]) # radii
# y_axis = []

radius = 5
ratio = 0.95
area = 0.75
baseline = np.sqrt(radius**2/(0.25+0.25*ratio**2))
# baseline = 2
# ratio = radius/0.577/baseline
temperature = 50
# ------------ iteratorrr
# mins = np.array([5, 6, 4, ])
# maxs = np.array([17, 18.5, 20, 21.5])
times = np.array([1,2,3,4,5,6,7])*86400
x_axis = times/66400 # radii
print('axis values', x_axis)
# y_axis = []
#--------------

for n, setting in enumerate(experiment_settings):
    scenario = setting[2]
    bus.data.options.set_manual(technique=setting[0])
    bus.data.options.set_manual(n_apertures=setting[1])
    bus.data.options.set_manual(array_config=setting[2])
    if figure_bool:
        fom = np.zeros(len(x_axis))
        fom_ter = np.zeros(len(x_axis))
        fom_hab = np.zeros(len(x_axis))
    if save_bool:
        path = root+'discussion'
        directory = os.path.join(path, scenario+'_'+what_are_we_varying)
        scenario_file = os.mkdir(directory)
    print('this run simulates an instrument which uses technique '+str(setting[0])+' with '+str(setting[1])+' apertures'
          ' in a(n) '+setting[2]+' config to loop over 500 close-by planets. Now lets create this instrument')
    for m, timer in enumerate(times):

        # baseline = np.sqrt(radius ** 2 / (0.25 + 0.25 * ratio ** 2))
        print('which yields for '+str(m)+' boundary values baseline: ', baseline, ', ratio: ', ratio,
              ', area: ', area, ' timer', timer)
        bus.data.options.set_manual(baseline=baseline,
                                    ratio=ratio,
                                    area=area,
                                    integration_time=timer)
        # create instrument
        instrument.apply_options()
        # fg.configuration_render(configuration=bus.data.inst['configuration'],
        #                         diameters=bus.data.inst['diameters'], save=figure_bool,
        #                         save_name='config')
        # --- check instrument
        mirror_radius = (area / bus.data.options.array['n_apertures'] / np.pi) ** 0.5
        dimension_check, _ = baseline_distances(bus.data.inst['configuration'])
        if mirror_radius >= dimension_check*0.5:
            print('---- rejected instrument')
            continue
        instrument.get_snr(verbose=False)
        # -------- Figure configuration
        # fg.configuration_render(bus.data.inst['configuration'], bus.data.inst['diameters'], save=False,
        #                         save_name=scenario+' configuration a_e_e')
        # wl_slab = 0
        # instrument.run_socket(s_name='transmission', method='plot_psf', slab=wl_slab, save=False,
        #                       save_name=scenario+' psf plot')
        # -------- Figure of all
        if figure_bool:
            fom[m] = len(bus.data.catalog[bus.data.catalog.snr_1h > 5])
            fom_ter[m] = len(bus.data.catalog[(bus.data.catalog.radius_p < 2) & (bus.data.catalog.snr_1h > 5)])
            fom_hab[m] = len(bus.data.catalog[(bus.data.catalog.habitable == True) & (bus.data.catalog.snr_1h > 5)])

        print('~done!~')
        if save_bool:
            test_name = 'disc_'+str(radius)
            output_path = os.path.join(directory, test_name)
            bus.data.export_catalog_txt(output_path=output_path+'.txt')
            bus.data.export_catalog(output_path=output_path+'.hdf5')
            bus.data.export_instrument(output_path=output_path)

# print('end of render @', time.perf_counter())


# ------------------  figure settings

figure_location = 'C:/Users/rogie/Desktop/Afstudeerproject/05. Figures (made by me)/02. Thesis/'


# ------------------ initiate figure
if figure_bool:
    title = experiment_settings[0][2]+what_are_we_varying+' sensitvity analysis '
    fig, ax1 = plt.subplots()
    ax1.plot(x_axis, fom, linewidth=2, marker='x', markersize=5, label='Total planet detections')
    ax1.plot(x_axis, fom_ter, linewidth=2, marker='x', markersize=5, label='Terrestrial detections')
    ax1.plot(x_axis, fom_hab, linewidth=2, marker='x', markersize=5, label='Habitable detections')
    # ax2.plot(baselines_emma, circle_fom, linewidth=2, marker='x', markersize=5, label='Total planet detections')
    # ax2.plot(baselines_emma, circle_fom_ter, linewidth=2, marker='x', markersize=5, label='Terrestrial detections')
    # ax2.plot(baselines_emma, circle_fom_hab, linewidth=2, marker='x', markersize=5, label='Habitable detections')
    ax1.legend()
    # ax2.legend()
    # ax1.set_title()
    ax1.set_ylabel('Detections [-]',fontsize =18)
    ax1.set_xlabel('Observation time [Days]',fontsize =18)
    ax1.grid()

    if save_bool:
        pic_name = title
        plt.savefig(figure_location+pic_name, dpi=500)
    plt.show()
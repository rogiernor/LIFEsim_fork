import lifesim
import numpy as np
# import matplotlib.patches as ptch
import matplotlib.pyplot as plt
from lifesim.core.modules import PhotonNoiseModule, TransmissionModule
# from matplotlib import colors
# import itertools
import lifesim.util.figures as fg
import os
# import sympy as sp
# import scipy as sc
import cmasher as cmr
# ------- CREATE INSTRUMENT
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
# connect all modules
bus.connect(('inst', 'transm'))
bus.connect(('inst', 'exo'))
bus.connect(('inst', 'local'))
bus.connect(('inst', 'star'))
bus.connect(('inst', 'ins_shot'))
bus.connect(('star', 'transm'))
# --------------- import test pop
root = 'C:/Users/rogie/Desktop/Afstudeerproject/03. Code/01. LIFEsim-master/'
planet_path = os.path.join(root, 'docs')
bus.data.import_catalog(
    input_path=planet_path+'/FOM_ppop.hdf5')
bus.data.catalog_remove_distance(stype=0, mode='larger', dist=0.)  # remove all A stars
bus.data.catalog_remove_distance(stype=4, mode='larger', dist=5.)  # remove M stars > 10pc
bus.data.catalog_remove_distance(stype=3, mode='larger', dist=5.)  # remove K stars > 10pc
bus.data.catalog_remove_distance(stype=2, mode='larger', dist=5.)  # remove G stars > 10pc
bus.data.catalog_remove_distance(stype=1, mode='larger', dist=5.)  # remove F stars > 10pc
# --------------- Select population of 50 planets
print(len(bus.data.catalog))
ppop_pop_size = 500
ppop_parser = len(bus.data.catalog)//ppop_pop_size
print('the amount of skips', ppop_parser)
bus.data.catalog = bus.data.catalog.loc[::ppop_parser]
# bus.data.catalog = bus.data.catalog.loc[bus.data.catalog['nuniverse'] == 0]
# ---- Set options
bus.data.options.set_scenario()
# ----- Set fr Options
max_rad = 4.0
x = [7.53894,	1.33138
, 2, 3]

# x[1] = max_rad/( 0.577 * x[0])
# x[1] = np.sqrt(max_rad**2/(0.25*x[0]**2)-1)
# print(x[1])
test = 'diamond erxtra defense'
save_bool = False
fig_save_bool = True
plot_maps = True
wl_slab = 18
# ------- Run with it
experiment_settings = [[2, 3, 'TTN'],
                       [1, 4, 'angel woolf'],
                       [1, 4, 'emma x'],
                       [1, 4, 'diamond'],
                       [2, 4, 'kite'],
                       [2, 5, 'pentagon']]
bus.data.options.set_manual(n_apertures=experiment_settings[x[3]][1],
                            technique=experiment_settings[x[3]][0],
                            array_config=experiment_settings[x[3]][2],
                            baseline=x[0],
                            ratio=x[1],
                            area=x[2],
                            instrument_temperature=50)

instrument.apply_options()
limit = 7
fg.configuration_render(configuration=bus.data.inst['configuration'],
                        diameters=bus.data.inst['diameters'], dimensions=limit,
                        add_circle=True,
                        save=fig_save_bool, save_name=test+' config')
maps = np.swapaxes(bus.data.inst['maps'], 0, 1)
extent_value = bus.data.inst['hfov_mas'][wl_slab]
extent = [-extent_value, extent_value, -extent_value, extent_value]
if plot_maps:
    # fg.ttn_figures(maps,
    #                 slab=wl_slab,
    #                 extent=extent,
    #                 save=save_bool,
    #                 save_name=test+ 'maps')
    # kernels = [maps[0], maps[1] - maps[2], maps[3] - maps[4], maps[5] - maps[6]]
    # # kernels = [maps[0], maps[1] - maps[2], maps[3] - maps[4], maps[5] - maps[6], maps[7] - maps[8], maps[9] - maps[10], maps[11] - maps[12]]
    # fg.transmission_figures(kernels,
    #                 slab=wl_slab,
    #                 extent=extent,
    #                 row_col=(2,2),
    #                 cmap=cmr.redshift,
    #                 save=save_bool,
    #                 save_name=test+'kernels')
    # fg.transmission_figures(maps,
    #                         slab=wl_slab,
    #                         row_col=(2, 2),
    #                         extent=extent,
    #                         save=fig_save_bool,
    #                         save_name=test+ 'maps')
    fg.single_map((maps[2])[wl_slab], extent=extent, labels=['alpha [mas]', 'beta[mas]'],
                  plot_title='$T_3$',
                  save=fig_save_bool, cmap='hot',
                  save_name=test + ' t_3', intensity='(norm.) transmission [-]')
    fg.single_map((maps[3])[wl_slab], extent=extent, labels=['alpha [mas]', 'beta[mas]'],
                  plot_title='$T_4$',
                  save=fig_save_bool, cmap='hot',
                  save_name=test + ' t_4', intensity='(norm.) transmission [-]')
    fg.single_map((maps[2]-maps[3])[wl_slab], extent=extent, labels = ['alpha [mas]', 'beta[mas]'], plot_title='Emma X-array chop',
                  save=fig_save_bool,
                  save_name=test+ ' chop', intensity='Modulation [-]')
instrument.run_socket(s_name='transmission', method='plot_psf', slab=wl_slab,
                      save=fig_save_bool, save_name=test+' rot psf', sliced=True, custom=False)
# ------ test FOM
instrument.get_snr(verbose=False)
fom_larger = bus.data.catalog[bus.data.catalog.snr_1h > 5.]
print(fom_larger[['radius_p', 'temp_p', 'distance_s', 'radius_s', 'habitable', 'snr_1h']]) #, 'dista
print(bus.data.catalog[bus.data.catalog.habitable == True].snr_1h.max())
print(bus.data.catalog[bus.data.catalog.habitable == True].snr_1h.idxmax())
print(len(fom_larger))
fom_ter = len(bus.data.catalog[(bus.data.catalog.radius_p < 2) & (bus.data.catalog.snr_1h > 5)])
fom_hab = len(bus.data.catalog[(bus.data.catalog.habitable == True) & (bus.data.catalog.snr_1h > 5)])
print('secondary obj,', fom_ter, fom_hab)
if save_bool == True:
    output_path='C:/Users/rogie/Desktop/Afstudeerproject/09. Experiments/02. Output/Champions!/'
    bus.data.export_catalog_txt(output_path=output_path+test+'.txt')
    bus.data.export_catalog(output_path=output_path + test + '.hdf5')
    bus.data.export_instrument(output_path=output_path+test)
# rho_list = []
# tips = np.array([10, 1000, 10000])*lifesim.util.constants.rad_per_mas
# for i, tip in enumerate(tips):
#     tilt = 0
#     focus = 100
#     hfov = 2 * np.arctan(bus.data.inst['diameter']/(2 * focus))        # np.array(bus.data.inst['hfov'])
#     print('hfov', hfov)
#     # if hfov.shape[-1] > 1:
#     #     print('RESHAPE!')
#     #     hfov = np.reshape(hfov, (hfov.shape[-1], 1, 1))
#     #     print('hfov shape', np.shape(hfov))
#     image_size = bus.data.options.other['image_size']
#     # recreate x and y but add tip and tilt angles
#     angle = np.linspace(-1, 1, image_size)
#     alpha = np.tile(angle, (image_size, 1))
#     beta = alpha.T
#     alpha = alpha * hfov
#     beta = beta * hfov
#     print('"pixel" size', alpha[0, 0] - alpha[0, 1])
#     rho, a_couple = instrument.run_socket(method='coupling_efficiency_calculator', s_name='transmission',
#                                           m=0, alpha=alpha, beta=beta, tip=tip, tilt=tilt, focus=focus, iterator=i)
#     rho_sum_square = np.sqrt(sum(np.array(rho)**2))
#     print(rho_sum_square)
#     rho_list.append(rho_sum_square)
# # coupling efficiency figure
# fig, ax = plt.subplots()
# ax.semilogx(tips, rho_list)
# ax.set_ylabel('coupling efficiency')
# ax.set_xlabel('tip [rad]')
# ax.set_ylim([0, 100])
# plt.show()

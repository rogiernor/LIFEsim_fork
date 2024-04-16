import lifesim
import matplotlib.patches as ptch
import matplotlib.pyplot as plt
from matplotlib import colors
import itertools
import numpy as np
import lifesim.util.figures as fg
import os
import sympy as sp
import scipy as sc
print(lifesim.util.constants.c)  # run test
import os

# ----------------------- Used to check validity of Optimization champions
# --- Modules
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
# --- connect all modules
bus.connect(('inst', 'transm'))
bus.connect(('inst', 'exo'))
bus.connect(('inst', 'local'))
bus.connect(('inst', 'star'))
bus.connect(('inst', 'ins_shot'))
bus.connect(('star', 'transm'))
# --- set scenario
bus.data.options.set_scenario()
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
# -------------- instrument
# ----- set manual
# bus.data.options.set_manual(baseline=2.0, ratio=3.2, area=4, integration_time=86400, instrument_temperature=37.19)
# ---- Classic settings
# geometry = 'emma x'
# randopm config
# x = [5.98311,5.59534,0.76875,5.99424,1.48339,2.49179,2.43732,4.20814]
# x1 = [x[0] * np.cos(x[4]), x[0] * np.sin(x[4])]
# x2 = [x[1] * np.cos(x[5]), x[1] * np.sin(x[5])]
# x3 = [x[2] * np.cos(x[6]), x[2] * np.sin(x[6])]
# x4 = [x[3] * np.cos(x[7]), x[3] * np.sin(x[7])]
# -------------- settings

loose = 0
# xes = [[6.27716,	1.4833,	2.49179]]
# xes = [[0.57296439, 5.64322003, 4.88063418],
#        [4.92170999, 2.59205801, 0.21470467],
#        [1.29566215, 5.86571193, 2.76487833]]
# xes = [[1.82972148, 4.39185359, 2.89754538, 0.19328092],
# [5.40652884 ,5.07771704 ,4.00080411, 1.68938942],
# [4.4985754,  5.04387415, 0.58308469 ,3.25564848]]
# x = [[0.892967,	4.92171,	2.59206,	0.214705,	3.9209],
# [1.62419,	3.43192,	2.55919,	1.11203,	6.09238],
# [4.16378,	5.47311,	5.09709,	2.35069,	0.47249]]
# xes = [[1, 2, 4.5, 5]]

xes = [[4.36571,	0.913388,	2.84737]]

# base_radii = [2.12899412, 2.82956014, 3.5766949]
# base_radii = [3.57363856, 2.32195104, 4.32510123]
# base_radii = [2.12899412, 2.82956014, 3.5766949]
base_radii = [5, 5, 5]
area = 2

tests = ['Random 3 ap winner 1']
save_bool = True
slice_bool = False
# ----------------- Test
for m, x in enumerate(xes):
    dim = len(x)
    test = tests[m]
    r = [base_radii[m]] * dim
    if loose >= 1:
        r[-loose:] = x[-loose:]
    geometry_list = []
    for i in list(range(dim)):
        aperture_position = [r[i] * np.cos(x[i]), r[i] * np.sin(x[i])]
        geometry_list.append(aperture_position)
    geometry = np.array(geometry_list)
    print(geometry)
    bus.data.options.set_manual(area=area, n_apertures=len(geometry), technique=2, array_config=geometry)
    bus.data.options.set_manual(wl_min=4., wl_max=18.5)
    bus.data.options.set_manual(spec_res=20)
    # ----- create instrument
    instrument.apply_options()
    # print('area', bus.data.inst['telescope_area'])
    # ------ figures of instrument
    wl_slab = 18
    extent_value = bus.data.inst['hfov_mas'][wl_slab]
    extent=[-extent_value, extent_value, -extent_value, extent_value]
    maps = np.swapaxes(bus.data.inst['maps'], 0, 1)
    # maps = bus.data.inst['norm_maps']
    # # print('min value maps', np.amin(bus.data.inst['norm_maps'][2, 0]))
    if len(geometry) == 3:
        fg.ttn_figures(maps, slab=wl_slab, extent=extent, save=save_bool, save_name=test+' TTN winner normalized')
    else:
        kernels = [maps[0], maps[1]-maps[2], maps[3]-maps[4], maps[5]-maps[6]]
        fg.transmission_figures(kernels,
                                slab=wl_slab,
                                row_col=(2, 2),
                                extent=extent,
                                save=save_bool,
                                save_name=test+'kernels')
    # # fg.ttn_figures(maps, slab=wl_slab, extent=extent)
    limit = max(r)
    fg.configuration_render(configuration=bus.data.inst['configuration'],
                            diameters=bus.data.inst['diameters'], dimensions=6, save=save_bool, save_name=test+' Winner config', add_circle=True)
    # # fg.plot_outputs_smart(bus.data.inst['U'])
    instrument.run_socket(s_name='transmission', method='plot_psf', slab=wl_slab,
                          save=slice_bool, save_name=test + ' rot psf', sliced=True, custom=False)
    instrument.get_snr(verbose=False)
    fom_larger = bus.data.catalog[bus.data.catalog.snr_1h > 5.]
    fom1 = len(bus.data.catalog[bus.data.catalog.snr_1h < 5])
    fom2 = len(bus.data.catalog[(bus.data.catalog.radius_p < 2) & (bus.data.catalog.snr_1h < 5)])
    print(fom_larger[['radius_p', 'sep_p', 'distance_s', 'temp_p', 'habitable', 'snr_1h']]) #, 'distance_s', 'temp_p', 'habitable', 'snr_1h'
    print('FOM len', len(fom_larger))
    print('FOM1 len', fom1)
    print('FOM2 len', fom2)
    if save_bool == True:
        output_path='C:/Users/rogie/Desktop/Afstudeerproject/09. Experiments/02. Output/Champions!/'
        bus.data.export_catalog_txt(output_path=output_path+test+'.txt')
        bus.data.export_catalog(output_path=output_path + test + '.hdf5')
        bus.data.export_instrument(output_path=output_path+test)

import lifesim
import time
import numpy as np
import lifesim.util.figures as fg
import os
import matplotlib.pyplot as plt
from lifesim.util.normalizer import baseline_distances
import lifesim.util.radiation as rd

# bus
bus = lifesim.Bus()
# modules
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
# options
bus.data.options.set_scenario()
root = 'C:/Users/rogie/Desktop/Afstudeerproject/09. Experiments/02. Output/'
figure_file = 'C:/Users/rogie/Desktop/Afstudeerproject/05. Figures (made by me)/02. Thesis/'
# -------------
# planet_path = "/docs/FOM_ppop.hdf5"
# bus.data.import_catalog(
#     input_path=planet_path)
# print('done',  len(bus.data.catalog))
# bus.data.catalog_remove_distance(stype=0, mode='larger', dist=0.)  # remove all A stars
# bus.data.catalog_remove_distance(stype=4, mode='larger', dist=5.)  # remove M stars > 10pc
# bus.data.catalog_remove_distance(stype=3, mode='larger', dist=5.)  # remove K stars > 10pc
# bus.data.catalog_remove_distance(stype=2, mode='larger', dist=5.)  # remove G stars > 10pc
# bus.data.catalog_remove_distance(stype=1, mode='larger', dist=5.)  # remove F stars > 10pc
# # --------------- Select population of 50 planets
# bus.data.catalog.info()
# ppop_pop_size = 1
# ppop_parser = len(bus.data.catalog)//ppop_pop_size
# # print('the amount of skips', ppop_parser)
# bus.data.catalog = bus.data.catalog.loc[::ppop_parser]
# print(bus.data.catalog[['radius_p', 'temp_p','distance_s','temp_s','z']])
# ----------------==  EArth clone
Earth_clone = np.array([1.0, 365, 1.0, 0.0, 0, 0, 0, 0., 0.29,
                        0, 0, 3.0, 1.0, 1.0, 0.4, 0, 0, 0,
                        285, 1.0, 0, 5778, 2.5, 0.0, 0., 0, 0, 0, 0]).reshape((1, 29))

bus.data.catalog_add_planet(Earth_clone)
bus.data.catalog.loc[0, 'lon'] = 135.           # this is the proper format for indexing dataframes
bus.data.catalog.loc[0, 'lat'] = 45.
print('what is sep', bus.data.catalog.loc[0, 'sep_p'])
bus.data.catalog.loc[0, 'angsep'] = bus.data.catalog.loc[0, 'sep_p']/bus.data.catalog.loc[0, 'distance_s']
print(bus.data.catalog.iloc[0])
print(bus.data.catalog.info)
# ------------------ Bruh
setting = [1, 4, 'emma x']
test_name = 'x-array temperature change large area check'
temperatures = [40, 50, 80]
save_bool = True
# --------------
print('start at', time.perf_counter())
scenario = setting[2]
bus.data.options.set_manual(technique=setting[0])
bus.data.options.set_manual(n_apertures=setting[1])
bus.data.options.set_manual(array_config=setting[2])
bus.data.options.set_manual(area=2, baseline=2, ratio=6)
fig, ax1 = plt.subplots()
# ax2 = ax1.twinx()
for t, temperature in enumerate(temperatures):
    print('which yields for iteration', t, 'temperature', temperature)
    bus.data.options.set_manual(instrument_temperature=temperature)
    # create instrument
    instrument.apply_options()
    # --- check instrument
  # earth clone: this array gets fed into data.single[]
    # ['03', 'h2o', 'co2']
    # print('shape of boys',  np.shape(bus.data.inst['wl_bins']),
    #                             np.shape(bus.data.inst['wl_bin_widths']),
    #                                      np.shape( bus.data.catalog['temp_p']),
    #                             np.shape(bus.data.catalog['radius_p']),
    #                                      np.shape( bus.data.catalog['distance_s']))
    planet_bb = rd.black_body(mode='planet',
                                bins=bus.data.inst['wl_bins'],
                                width=bus.data.inst['wl_bin_widths'],
                                temp=bus.data.catalog['temp_p'].values,
                                radius=bus.data.catalog['radius_p'].values,
                                distance=bus.data.catalog['distance_s'].values)
    #
    #
    snr_specs, signals, noises, noise_bg_list = instrument.get_spectrum(temp_s=bus.data.catalog['temp_s'].values,
                                                                        radius_s=bus.data.catalog['radius_s'].values,
                                                                        distance_s=bus.data.catalog['distance_s'].values,
                                                                        lat_s=bus.data.catalog['lat'].values,
                                                                        z=bus.data.catalog['z'].values,
                                                                        angsep=bus.data.catalog['angsep'].values,
                                                                        input_flux_planet_spectrum=planet_bb,
                                                                        integration_time=bus.data.options.other
                                                                        ['integration_time'],
                                                                        lines_bool=True,
                                                                        gasses=['o3'])
    # print('signall', signals)
    # wl_slab=21
    # print('weird slab', bus.data.inst['wl_bins'][wl_slab])
    # instrument.run_socket(s_name='transmission', method='plot_psf', slab=wl_slab,
    #                       save=save_bool, save_name='emma x rot psf', sliced=True, custom=True)
    # instrument.get_snr(False)

    ax1.semilogy(bus.data.inst['wl_bins']*1e6, bus.data.inst['ti_leak'], linestyle='--', label = '$T_i$ ='+str(temperature)+' K')
    print('~done!~')
    # print('amount of bins', len(snr_specs[0]))
    # if save_bool:
    #     test_name = '557_5_pc_planets_'+str(baseline)+'_'+str(ratio)
    #     output_path = os.path.join(directory, test_name)
    #     bus.data.export_catalog_txt(output_path=output_path+'.txt')
    #     bus.data.export_catalog(output_path=output_path+'.hdf5')
    #     bus.data.export_instrument(output_path=output_path)

ax1.semilogy(snr_specs[0]*1e6, noise_bg_list[0], label='ez leak')
ax1.semilogy(snr_specs[0]*1e6, noise_bg_list[1], label = 'lz noise')
ax1.semilogy(snr_specs[0]*1e6,  noise_bg_list[2], label = 'sl noise')
# ax1.set_ylabel('Noise flux [ph/s/bin]', fontsize=18)
ax1.set_ylabel('Noise [ph/s/bin]', fontsize=18)
ax1.set_xlabel('$\lambda$ [$\mu$m]', fontsize=18)
# ax1.set_ysa([])
ax1.set_ylim([1e-1, 1e5])
ax1.grid()
# ax1.legend(loc='upper left')
ax1.legend(loc='lower right', fontsize=14)
fig.tight_layout()
if save_bool:
    plt.savefig(figure_file+test_name+'.png')
plt.show()
# print('end of render @', time.perf_counter())

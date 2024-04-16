import lifesim
import os
# import matplotlib.patches as ptch
import matplotlib.pyplot as plt
import numpy as np
import lifesim.util.figures as figures
from tqdm import tqdm
from lifesim.util.normalizer import baseline_distances
print(lifesim.util.constants.c)


# ---------- Set-Up ----------

# create bus
bus = lifesim.Bus()
# obligate simulation settings
# ---------- Loading the Catalog ----------
root = 'C:/Users/rogie/Desktop/Afstudeerproject/03. Code/01. LIFEsim-master/'
bus.data.import_catalog(
    input_path= root+'/docs/FOM_ppop.hdf5')
print(''
      ''
      'done'
      ''
      ''
      '')
bus.data.catalog_remove_distance(stype=0, mode='larger', dist=0.)  # remove all A stars
bus.data.catalog_remove_distance(stype=4, mode='larger', dist=5.)  # remove M stars > 10pc to
bus.data.catalog_remove_distance(stype=3, mode='larger', dist=5.)  # remove M stars > 10pc to
bus.data.catalog_remove_distance(stype=2, mode='larger', dist=5.)  # remove M stars > 10pc to
bus.data.catalog_remove_distance(stype=1, mode='larger', dist=5.)  # remove M stars > 10pc to
# speed up calculation -> make selection
bus.data.catalog.info()
ppop_pop_size = 500
ppop_parser = len(bus.data.catalog)//ppop_pop_size
# print('the amount of skips', ppop_parser)
bus.data.catalog = bus.data.catalog.loc[::ppop_parser]
# ---------- Creating the Instrument ----------

# create modules and add to bus
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

# connect all modules
bus.connect(('inst', 'transm'))
bus.connect(('inst', 'exo'))
bus.connect(('inst', 'local'))
bus.connect(('inst', 'star'))
bus.connect(('star', 'transm'))

# ----- test settings
# parent_name = 'design space output'
# test_name = scenario+'_baseline_test'  # TODO this is important, brother
# new_dir = os.path.join(root, parent_name, test_name)
# os.mkdir(new_dir)

# ----------- monte charlo test
n_instances = 100
save_bool = True
np.random.seed(bus.data.options.other['random_seed'])
for i in tqdm(range(100)):
    # Syst settings
    area = 1
    x = np.random.random(4)*2*np.pi
    r = np.random.random(1)*6
    print('íteration', i)
    print('azimuthal coordinates', x)
    print('radius', r)
    x1 = [r * np.cos(x[0]), r * np.sin(x[0])]
    x2 = [r * np.cos(x[1]), r * np.sin(x[1])]
    x3 = [r * np.cos(x[2]), r * np.sin(x[2])]
    x4 = [r * np.cos(x[3]), r * np.sin(x[3])]
    geometry = np.array([x1, x2, x3, x4])
    bus.data.options.set_manual(n_apertures=4, technique=2, array_config=geometry)
    instrument.apply_options()
    # check dimensions

    mirror_radius = (area / bus.data.options.array['n_apertures'] / np.pi) ** 0.5
    dimension_check, _ = baseline_distances(bus.data.inst['configuration'])
    if mirror_radius >= dimension_check * 0.5:
        print('---- rejected instrument')
        continue
    # -- get_snr
    instrument.get_snr(verbose=False)
    # # ----- show fom
    fom_larger = bus.data.catalog[bus.data.catalog.snr_1h > 5.]
    print(fom_larger[['radius_p', 'sep_p', 'distance_s', 'temp_p', 'habitable', 'snr_1h']])
    print(fom_larger.snr_1h.values)
    print(fom_larger.sep_p.values)
    # # ------ Save resulting
    # output_path = os.path.join(new_dir, str(instance))
    if save_bool:
        path = 'C:/Users/rogie/Desktop/Afstudeerproject/09. Experiments/Experiment data/kernel nulling MC'
        run_name = 'mc2'
        directory = os.path.join(path, run_name)
        scenario_file = os.mkdir(directory)
        instance_name = 'random_config_'+str(i)
        output_path = os.path.join(directory, instance_name)
        bus.data.export_catalog(output_path=output_path+'.hdf5')
        bus.data.export_catalog_txt(output_path=output_path + '.txt')
        bus.data.export_instrument(output_path=output_path + '_instrument')
print('simulation finished')



from warnings import warn

import numpy as np
import sympy as sp
from pathlib import Path
from tqdm import tqdm
from spectres import spectres
from PyQt5.QtGui import QGuiApplication

from lifesim.core.modules import InstrumentModule
from lifesim.util.habitable import single_habitable_zone
from lifesim.util.radiation import black_body, spectral_lines
from lifesim.util.normalizer import transmission_normalizer, onedimensional_map_normalizer
import lifesim.util.figures as fg
import matplotlib.pyplot as plt


class Instrument(InstrumentModule):
    """
    The Instrument class represents the central module for simulating the LIFE array. It connects
    to other modules which calculate signal and noise terms and distributes tasks and data
    between them. The instrument class features two socket types:
        a)  For calculation of the instrument transmission map a single socket of f_type
            'transmission'.
        b)  For simulation of the photon noise sources a number (set in the options class) of
            sockets of f_type 'photon_noise'.

    Notes
    -----
    Note, that all attributes are saved in the data class.

    Attributes
    ----------
    data.options : lifesim.Options
        The options class containing all setting for the array and computations.
    data.inst['bl'] : float
        Length of the shorter, nulling baseline in [m].
    data.inst['telescope_area'] : float
        Area of all array apertures combined in [m^2].
    data.inst['eff_tot'] : float
        Total efficiency of the telescope as ratio of generated counts over incoming photons
        (dimensionless).
    data.inst['wl_bins'] : np.ndarray
        Central values of the spectral bins in the wavelength regime in [m].
    data.inst['wl_bin_widths'] : np.ndarray
        Widths of the spectral wavelength bins in [m].
    data.inst['wl_bin_edges'] : np.ndarray
        Edges of the spectral wavelength bins in [m]. For N bins, this array will contain N+1
        edges.
    data.inst['hfov'] : np.ndarray
        Contains the half field of view of the observatory in [rad] for each of the spectral bins.
    data.inst['hfov_mas'] : np.ndarray
        Contains the half field of view of the observatory in [milliarcseconds] for each of the
        spectral bins.
    data.inst['rad_pix'] : np.ndarray
        Contains the size of each pixel projected to the sky in [rad].
    data.inst['mas_pix'] : np.ndarray
        Contains the size of each pixel projected to the sky in [milliarcseconds].
    data.inst['apertures'] : np.ndarray
        Positions of the collector spacecraft relative to the beam combiner in [m].
    data.inst['radius_map'] : np.ndarray
        A map used for speeding up calculations. Contains the distance of a pixel
        from the center of the detector in [pix].
    """

    def __init__(self,
                 name: str):
        """
        Parameters
        ----------
        name : str
            Name of the instrument module.
        """

        super().__init__(name=name)

    def apply_options(self, verbose: bool = False, zoom: float = 1.):
        """
        Applies the options given to the instrument module and recalculates all necessary values.
        it generates an instrument, including the maps.

        :param
            none, formally but,
            verbose: bool = False
                If true, apply_options will spit out its options.
                Which at the moment is only a figure of the first 4 maps
            zoom: float = 1.
                Method of decoupling the FOV from the size of the graph

        :returns
            None
        """
        #   #    Important
        # --------- Import from Options & called up top TODO SETTINGS -> MOVE BACK TO OPTIONS?
        self.data.inst['bl'] = self.data.options.array['baseline']
        self.data.inst['n_apertures'] = self.data.options.array['n_apertures']
        self.data.inst['array_config'] = self.data.options.array['array_config']
        self.data.inst['quantum_eff'] = self.data.options.array['quantum_eff']
        self.data.inst['throughput'] = self.data.options.array['throughput']
        self.data.inst['combination_technique'] = self.data.options.array['technique']
        self.data.inst['integration_time'] = self.data.options.other['integration_time']
        self.data.inst['fiber_angle'] = self.data.options.other['fiber_angle']
        # - ------------------------- Get array parameters from options for faster calculation
        self.data.inst['diameter'] = np.sqrt(self.data.options.array['area'] / np.pi / self.data.inst['n_apertures'])\
                                     * 2.
        self.data.inst['diameters'] = np.repeat(self.data.inst['diameter'], self.data.inst['n_apertures'])
        self.data.inst['telescope_area'] = self.data.options.array['area']
        self.data.inst['eff_tot'] = self.data.inst['quantum_eff'] \
                                    * self.data.inst['throughput']
        # ---- Array configuration
        if isinstance(self.data.inst['array_config'], str):
            self.preset_configuration(case=self.data.inst['array_config'])
        else:
            self.data.inst['configuration'] = self.data.inst['array_config'] # TODO ADD SEED, SO SETTINGS CAN RE-CREATE
        # ----- Beam combination
        if self.data.inst['combination_technique'] == 2:
            self.kernel_beam_combiner(self.data.inst['n_apertures'])
        else:
            self.data.inst['U'] = self.beam_combination_setting(self.data.inst['array_config'])
        # -   - ----------------  Calculate the spectral channels with a constant spectral resolution
        self.data.inst['wl_bins'], \
        self.data.inst['wl_bin_widths'], \
        self.data.inst['wl_bin_edges'] = self.get_wl_bins_const_spec_res()

        # TODO remove the double usage of mas and rad, stick to only one

        # ---------------------  Pixel and FOV size
        self.data.inst['zoom'] = zoom
        self.data.inst['hfov'] = self.data.inst['zoom'] * self.data.inst['wl_bins'] \
                                 / (2. * self.data.inst['diameter'])

        self.data.inst['hfov_mas'] = self.data.inst['hfov'] * (3600000. * 180.) / np.pi
        self.data.inst['rad_pix'] = (2 * self.data.inst['hfov']) \
                                    / self.data.options.other['image_size']  # Radians per pixel
        self.data.inst['mas_pix'] = (2 * self.data.inst['hfov_mas']) \
                                    / self.data.options.other['image_size']  # mas per pixel

        # --------- ------------- coordinate maps for faster calculations
        x_map = np.tile(np.array(range(0, self.data.options.other['image_size'])),
                        (self.data.options.other['image_size'], 1))
        y_map = x_map.T
        r_square_map = ((x_map - (self.data.options.other['image_size'] - 1) / 2) ** 2
                        + (y_map - (self.data.options.other['image_size'] - 1) / 2) ** 2)
        self.data.inst['radius_map'] = np.sqrt(r_square_map)

        # ------------ Normalization settings?  ----------------  TODO here?
        self.data.inst['normalization_multiplier'] = 1

        # ------------ Generate transmission ------------------
        maps = self.run_socket(s_name='transmission',
                               method='transmission_map')

        # Non-normalized version for transmission_curve2 & pn star NB that 'norm' refers to non-normalized

        self.data.inst['norm_maps'] = maps
        multi = self.data.inst['normalization_multiplier']
        # normalize
        if self.data.inst['combination_technique'] == 2:
            self.data.inst['maps'] = multi * transmission_normalizer(transmission_map=maps,
                                                                     normalization_map=maps[0])

        elif self.data.inst['combination_technique'] == 1:
            self.data.inst['maps'] = multi * transmission_normalizer(transmission_map=maps,
                                                                     normalization_map=maps[0]*2)
            if verbose:
                scale_value = self.data.inst['hfov_mas'][11]
                fg.transmission_figures(maps=list(np.swapaxes(self.data.inst['maps'], 0, 1)), slab=11, row_col=(2, 2),
                                        extent=[-scale_value, scale_value,
                                                -scale_value, scale_value])
        else:
            self.data.inst['maps'] = multi * transmission_normalizer(transmission_map=maps,
                                                                     normalization_map=maps[0])
        # ------------- Generate fiber #model: ruilier 1999
        # numerical_aperture_fiber = np.sin(self.data.inst['fiber_angle'] * np.pi / 180)
        # fiber_diameter = 2.40499*self.data.inst['wl_bins'] / (2 * np.pi * numerical_aperture_fiber) # ORR self defined? longeur d'onde de cioupure
        # normalized_frequency = 2 * np.pi * fiber_diameter * numerical_aperture_fiber / self.data.inst['wl_bins']
        # self.data.inst['fundamental_mode_width'] = fiber_diameter * (0.65 + 1.619 / normalized_frequency ** (3/2) + 2.879 / normalized_frequency**6)
        # ------------- Create instrumental noise
        if self.data.options.models['instrumental_photon_noise']:
            self.data.inst['ti_leak'] = self.run_socket(s_name='instrument_photon_noise', method='inst_t_noise')
            # print('ti leak retrieved', self.data.inst['ti_leak'])

    def get_wl_bins_const_spec_res(self):
        """
        Create the wavelength bins for the given spectral resolution and wavelength limits.
        """
        wl_edge = self.data.options.array['wl_min']
        wl_bins = []
        wl_bin_widths = []
        wl_bin_edges = [wl_edge]

        while wl_edge < self.data.options.array['wl_max']:

            # set the wavelength bin width according to the spectral resolution
            wl_bin_width = wl_edge / self.data.options.array['spec_res'] / \
                           (1 - 1 / self.data.options.array['spec_res'] / 2)

            # make the last bin shorter when it hits the wavelength limit
            if wl_edge + wl_bin_width > self.data.options.array['wl_max']:
                wl_bin_width = self.data.options.array['wl_max'] - wl_edge

            # calculate the center and edges of the bins
            wl_center = wl_edge + wl_bin_width / 2
            wl_edge += wl_bin_width
            wl_bins.append(wl_center)
            wl_bin_widths.append(wl_bin_width)
            wl_bin_edges.append(wl_edge)

        # convert everything to [m]
        wl_bins = np.array(wl_bins) * 1e-6  # in m
        wl_bin_widths = np.array(wl_bin_widths) * 1e-6  # in m
        wl_bin_edges = np.array(wl_bin_edges) * 1e-6  # in m

        return wl_bins, wl_bin_widths, wl_bin_edges

    def get_snr(self,
                verbose: bool = False,
                var_file: int = None):
        """
        Calculates the signal-to-noise ration for all planets in the catalog. Also works for singular planets, although
        you might have to trick the method to think its an array of planets.

        Parameters
        ----------
        verbose : bool
            if verbose = True, do everything VERY LOUDLY
        """

        # options are applied before the simulation run
        if self.data.inst:
            pass
        else:
            raise ValueError('ERROR: no options applied before get_SNR')

        integration_time = self.data.inst['integration_time']
        self.data.inst['modulation_efficiency'] = []
        self.data.catalog['snr_1h'] = np.zeros_like(self.data.catalog.nstar, dtype=float)

        # create mask returning only unique stars (so planets aren't run ^2 times)
        _, temp = np.unique(self.data.catalog.nstar, return_index=True)

        star_mask = np.zeros_like(self.data.catalog.nstar, dtype=bool)
        star_mask[temp] = True

        # iterate over all stars
        for i, n in enumerate(tqdm(np.where(star_mask)[0])):
            nstar = self.data.catalog.nstar.iloc[n]

            # calculate the noise from the background sources
            noise_bg_list = self.run_socket(s_name='photon_noise',
                                            method='noise',
                                            index=n)
            # noise_bg_list.pop(0)
            # print('noise', noise_bg_list)
            # TODO: Reinstate the method in which the noise list is keyed by the name of the
            #  producing noise module
            # add instrumental photon noise
            if self.data.options.models['instrumental_photon_noise']:
                noise_bg_list.append(self.data.inst['ti_leak'])
            # print('noise_list', noise_bg_list)
            # sum up all noise
            if type(noise_bg_list) == list:
                noise_bg = np.zeros_like(noise_bg_list[0])
                for _, noise in enumerate(noise_bg_list):
                    noise_bg += noise
            else:
                noise_bg = noise_bg_list
            # print('noise before mult', noise_bg)

            noise_bg = noise_bg * integration_time * self.data.inst['eff_tot']
            # print('noiseafter mult', noise_bg)
            # go through all planets for the chosen star
            for _, n_p in enumerate(np.argwhere(
                    self.data.catalog.nstar.to_numpy() == nstar)[:, 0]):
                # print('n_p', n_p)
                # re-add exo_zodi
                # exo_zodi_noise = self.run_socket(s_name='photon_noise',
                #                                 method='noise2',
                #                                 index=n_p)
                #
                # # print('exo zodi noise', exo_zodi_noise[0])
                # exo_zodi_bg = exo_zodi_noise[0] * integration_time * self.data.inst['eff_tot']
                # print('noise after addition',  exo_zodi_bg)
                # calculate the photon flux originating from the planet
                flux_planet_thermal = black_body(mode='planet',
                                                 bins=self.data.inst['wl_bins'],
                                                 width=self.data.inst['wl_bin_widths'],
                                                 temp=self.data.catalog['temp_p'].iloc[n_p],
                                                 radius=self.data.catalog['radius_p'].iloc[n_p],
                                                 distance=self.data.catalog['distance_s'].iloc[n_p]
                                                 )
                # import matplotlib.pyplot as plt
                # plt.figure()
                # plt.plot(self.data.inst['wl_bins'], flux_planet_thermal)
                # plt.ylabel('Planet flux [$ph s^{-1} m^{-2} \mu m ^{-1}$')
                # plt.xlabel('$\lambda [\mu m]$')
                # plt.show()

                # calculate the transmission efficiency of the planets separation

                transm_eff, transm_noise = self.run_socket(s_name='transmission',
                                                           method='transmission_efficiency',
                                                           index=n_p)
                #print('ttansm_eff', transm_eff, 'transm_noise', transm_noise)
                # ----------- calculate the signal and photon noise flux received from the planet
                if self.data.inst['combination_technique'] == 2:
                    flux_planet = np.zeros((transm_eff.shape[0], transm_eff.shape[1]))
                    noise_planet = np.zeros((transm_eff.shape[0], transm_eff.shape[1]))
                    for m in range(transm_eff.shape[1]):
                        flux_planet_iterable = flux_planet_thermal \
                                               * transm_eff[:, m] \
                                               * integration_time \
                                               * self.data.inst['eff_tot'] \
                                               * self.data.inst['telescope_area']  # TODO DO NOT ADD WAVELENGTH BINS
                        flux_planet[:, m] = flux_planet_iterable
                        noise_planet_iterable = flux_planet_thermal \
                                                * transm_noise[:, m] \
                                                * integration_time \
                                                * self.data.inst['eff_tot'] \
                                                * self.data.inst['telescope_area']
                        noise_planet[:, m] = noise_planet_iterable
                else:
                    # print('not kernel')
                    # print('transmission efficiency', transm_eff)
                    # print('noise_efficiency', transm_noise)
                    # print('telescope area', self.data.inst['telescope_area'], self.data.inst['eff_tot'],
                    #       integration_time)
                    flux_planet = flux_planet_thermal \
                                  * transm_eff \
                                  * integration_time \
                                  * self.data.inst['eff_tot'] \
                                  * self.data.inst['telescope_area']
                    noise_planet = flux_planet_thermal \
                                   * transm_noise \
                                   * integration_time \
                                   * self.data.inst['eff_tot'] \
                                   * self.data.inst['telescope_area']

                # Add up the noise and calculate the SNR
                # print('constituent noises:', noise_bg, noise_planet, exo_zodi_bg)
                noise = (noise_bg + noise_planet) * 2
                # for bracewell it's only one noise per output
                if self.data.inst['combination_technique'] == 0:
                    noise = noise_bg + noise_planet
                # elif self.data.inst['combination_technique'] == 2:
                #     snr = flux_planet ** 2 / noise
                #     snr2 = np.swapaxes(snr, 0, 1)
                #     for i in snr2:
                #         snr_channel = np.sqrt(i.sum())
                #         print('snr per channel rsm', snr_channel)
                # print('check planets', n_p)


                self.data.catalog.snr_1h.iat[n_p] = np.sqrt(
                    (flux_planet ** 2 / noise).sum())
                if verbose:
                    if var_file:
                        print('transmission efficiency', transm_eff)
                        print('transmission noise', transm_noise)
                        print('noise list 1. ez_leak, 2. lz_noise, 3. sl_leak 4. instrumental', noise_bg_list[0][var_file], noise_bg_list[1][var_file],
                              noise_bg_list[2][var_file], noise_bg_list[3][var_file])
                        print('background noise', noise_bg[var_file])
                        print('planet noise', noise_planet[var_file])
                        print('total noise', noise[var_file])
                        print('signal planet', flux_planet[var_file])
                        print('SNR slab', np.sqrt(flux_planet[var_file]**2/noise[var_file]))
                        print('SNR for all wavelengths', np.sqrt(flux_planet ** 2 / noise))
                        print()
                        print('#')
                        print('Total SNR', np.sqrt((flux_planet ** 2 / noise).sum()))
                        print('#')
                        print()
                        save_bool = False
                        if self.data.inst['combination_technique'] == 2 and self.data.options.array['array_config'] != 'TTN' :
                            for blork in range(3):
                                print(blork)
                                labels = ['exozodi', 'localzodi', 'stellar', 'ínstrument', 'snr_kernel'+str(blork)]
                                snr = np.sqrt(flux_planet ** 2 / noise)
                                print('snr shape', np.shape(snr))
                                fg.plot_noise_snr( self.data.inst['wl_bins'], labels=labels, noises=np.moveaxis(noise_bg_list, -1, 0).tolist()[blork],
                                                  snr=np.swapaxes(snr, 0, 1).tolist()[blork], save=save_bool, save_name='kernel kerve'+str(blork))  # , snr=snr)
                        else:
                            labels = ['exozodi', 'localzodi', 'stellar', 'thermal', 'SNR', 'signal']
                            snr = np.sqrt(flux_planet ** 2 / noise)
                            fg.plot_noise_snr(self.data.inst['wl_bins'],  labels=labels,
                                              noises= noise_bg_list,
                                              snr=snr,
                                              signals=[flux_planet],
                                              save=save_bool, save_name='single output signal curve')
                    else:
                        import matplotlib.pyplot as plt
                        print('noise list 1. ez_leak, 2. lz_noise, 3. sl_leak 4. ', noise_bg_list[0][18], noise_bg_list[1][18],
                              noise_bg_list[2][18])
                        # print('flux signal planet', flux_planet)
                        # print('spectrum', flux_planet_thermal)
                        # print('noise total', noise_bg[0])
                        # print('noise intensity', noise_bg[0:5])
                        print('transmission efficiency', transm_eff[18])
                        print('transmission noise', transm_noise[18])
                        # print('thermal flux planet', sum(flux_planet_thermal**2), flux_planet_thermal.shape)
                        # print('efficiency over noise 0-10', transm_eff[0:6]/transm_noise[0:6])
                        # print('n star', nstar, 'n_planet', n_p)
                        # print('noise_planet', noise_planet[0], noise_planet.shape)
                        print('total noise', noise[18])
                        print('flux_planet', flux_planet[18])
                        # if nstar == 6:
                        save_bool = True
                        if self.data.inst['combination_technique'] == 2 and self.data.options.array['array_config'] != 'TTN' :
                            # print('nothing'
                            # )
                            labels = ['exozodi', 'localzodi', 'stellar', 'snr']
                            for blork in range(3):
                                print(blork)
                                # print('noise shape', noise_bg_list[blork])
                                labels = ['exozodi', 'localzodi', 'stellar', 'ínstrument', 'snr_kernel'+str(blork)]
                                snr = np.sqrt(flux_planet ** 2 / noise)
                                print('snr', np.shape(snr))

                                fg.plot_noise_snr(np.moveaxis(noise_bg_list, -1, 0).tolist()[blork], self.data.inst['wl_bins'], labels=labels,
                                                  snr=np.swapaxes(snr, 0, 1).tolist()[blork], save=save_bool, save_name='kernel kerve'+str(blork))  # , snr=snr)
                                # print('flux list', flux_planet)
                        else:
                            labels = ['exozodi', 'localzodi', 'stellar', 'thermal', 'SNR', 'signal', 'planet noise']
                            snr = np.sqrt(flux_planet ** 2 / noise)
                            # print('noise', noise_bg_list)m
                            fg.plot_noise_snr( self.data.inst['wl_bins'], labels=labels,
                                              noises=noise_bg_list,
                                              snr=snr,
                                              signals=[flux_planet, noise_planet],
                                              save=save_bool, save_name='all together')
                            labels = ['exozodi', 'localzodi', 'stellar', 'thermal', 'SNR', 'signal']
                            fg.plot_noise_snr(self.data.inst['wl_bins'], labels=labels,
                                              noises=noise_bg_list,
                                              save=save_bool, save_name='noise only')
                            labels = ['SNR', 'signal', 'planet noise']
                            print('type of sifgnals', type(flux_planet))
                            fg.plot_noise_snr(self.data.inst['wl_bins'], labels=labels,
                                              signals=[flux_planet, noise_planet],
                                              save=save_bool, save_name='signal only')

                            # LIFESIM Figure recreation
                            # labels = ['exozodi', 'localzodi', 'stellar', 'exozodi noise', 'localzodi noise', 'stellar noise', 'signal']
                            # print(np.shape(noise_bg_list[:3]))
                            # actual_noise = np.sqrt(noise_bg_list[:3])
                            # supernoise = np.vstack((noise_bg_list[:3], actual_noise))
                            # super_noise_sig = np.vstack((supernoise, flux_planet[np.newaxis, :]))
                            # print(np.shape(super_noise_sig))
                            # fg.plot_noise_snr(super_noise_sig, self.data.inst['wl_bins'], labels=labels)#, snr=snr)
                            # # fg.plot_noise_snr(noise_bg_list, self.data.inst['wl_bins'], labels=labels,
                            # #                   snr=flux_planet)#, snr=snr)
                            # print('flux list', flux_planet)
                            # fg.plot_noise_snr([flux_planet], self.data.inst['wl_bins'], labels=['signal', 'snr'], snr=snr)
                            # fig, ax1 = plt.subplots()
                            # ax1.stairs(snr, self.data.inst['wl_bin_edges']*10**6, color='k', label='snr')
                            # ax1.set_yscale('log')
                            # ax1.set_ylabel('SNR per bin [-]')
                            # ax1.legend(loc='upper right')
                            # ax2 = ax1.twinx()
                            # ax2.stairs(flux_planet, self.data.inst['wl_bin_edges']*10**6, linestyle='--', color='r', label='Planet')
                            # ax2.plot(self.data.inst['wl_bins']*10**6, np.sqrt(noise), linestyle='-.', color='b', label='Shot noise')
                            # ax2.set_ylabel('Detected signal per bin [e/s/bin]')
                            # ax2.set_xlabel('$\lambda$[$\mu$m]')
                            # ax2.set_yscale('log')
                            # plt.legend()
                            # plt.show()
                        print()
                        print('#')
                        print('SNR', np.sqrt((flux_planet ** 2 / noise).sum()))
                        print('#')
                        print()

    # TODO: fix units in documentation
    def get_spectrum(self,
                     temp_s: float,  # in K
                     radius_s: float,  # in R_sun
                     distance_s: float,  # in pc
                     lat_s: float,  # in radians
                     z: float,  # in zodis
                     angsep: float,  # in arcsec
                     input_flux_planet_spectrum: list,  # in ph m-3 s-1 over m
                     integration_time: float,  # in s
                     pbar=None,
                     baseline_to_planet: bool = False,
                     baseline: float = None,
                     safe_mode: bool = False,
                     lines_bool: bool = False,
                     gasses = None):
        """
        Calculate the signal-to-noise ratio per spectral bin of a given spectrum of a single
        planet.

        Parameters
        ----------
        temp_s : float
            Temperature of the observed star in [K].
        radius_s : float
            Radius of the observed star in [sun radii].
        distance_s : float
            Distance between the observed star and the LIFE array in [pc].
        lat_s : float
            Ecliptic latitude of the observed star in [rad].
        z : float
            Zodi level in the observed system in [zodis], i.e. the dust surface density of the
            observed system is z-times as high as in the solar system.
        angsep : float
            Angular separation between the observed star and the observed exoplanet in [arcsec].
        flux_planet_spectrum : list
            Spectrum of the planet. In the first element of the list `flux_planet_spectrum[0]`, the
            wavelength bins of the spectrum must be given in [m]. In the second element
            `flux_planet_spectrum[1]`, the photon count within the spectral bin must be given in
            [photons m-3 s-1].
        integration_time : float
            Time that the LIFE array spends for integrating on the observed planet in [s].
        pbar
            Takes a PyQt5 QProgressBar to display the progress of the baseline optimization.
        baseline_to_planet : bool
            If set to True, the baseline will be optimized to the position of the planet. If set
            to False, the baseline will be optimized to the center of the habitable zone of the
            host star.
        baseline : float
            Specifies a custom baseline. To have an effect, baseline_to_planet must be set to
            False.

        Returns
        -------
        Tuple[wl_bins, snr_spec]
            Returns the wavelength bins in [m] in the first element and the SNR per wavelength bin
            in the second element.
        flux_planet
            Returns the flux of the planet as used in the simulation in [photons]
        noise
            Returns the noise contribution in [photons]
        """

        # options are applied before the simulation run
        # self.apply_options()

        # write the given parameters to the single planet data in the bus. If the connected modules
        # are given an empty index to specify the star, they will use the data saved in this single
        # planet location
        if gasses is None:
            gasses = ['o3', 'h20', 'co2']
        self.data.single['temp_s'] = temp_s
        self.data.single['radius_s'] = radius_s
        self.data.single['distance_s'] = distance_s
        self.data.single['lat'] = lat_s
        self.data.single['z'] = z
        self.data.single['angsep'] = angsep

        # calculate the habitable zone of the specified star
        s_in, s_out, l_sun, \
        hz_in, hz_out, \
        hz_center = single_habitable_zone(model=self.data.options.models['habitable'],
                                          temp_s=temp_s,
                                          radius_s=radius_s)

        self.data.single['l_sun'] = l_sun

        # use spectres to rescale the spectrum onto the correct wl bins
        # flux_planet_spectrum_input = flux_planet_spectrum
        # print('input spectrum', flux_planet_spectrum_input)
        # print('new wavs', self.data.inst['wl_bin_edges'])
        # flux_planet_spectrum = spectres(new_wavs=self.data.inst['wl_bin_edges'],
        #                                 spec_wavs=flux_planet_spectrum[0],
        #                                 spec_fluxes=flux_planet_spectrum[1],
        #                                 edge_mode=True)

        if lines_bool:
            print('lines activated')
            line = 1
            for gas in gasses:
                line *= spectral_lines(self.data.inst['wl_bins'], gas=gas, thickness=-1)
            flux_planet_spectrum = input_flux_planet_spectrum * line
        else:
            flux_planet_spectrum = input_flux_planet_spectrum
        # print(flux_planet_spectrum, np.shape(flux_planet_spectrum))
        # plt.plot(self.data.inst['wl_bins'], line, label='lines_only')
        # plt.plot(self.data.inst['wl_bins'], flux_planet_spectrum, label='Effect of lines')
        # plt.plot(self.data.inst['wl_bins'], input_flux_planet_spectrum, linestyle='--', color='k', label='without gasses')
        # plt.legend()
        # plt.ylabel('Planet flux [$ph s^{-1} m^{-2} \mu m ^{-1}$')
        # plt.xlabel('$\lambda [\mu m]$')
        # plt.show()

        transm_eff, transm_noise = self.run_socket(s_name='transmission',
                                                   method='transmission_efficiency',
                                                   index=None)
        # print('tansm_eff', transm_eff, 'transm_noise', transm_noise)
        # calculate the signal and photon noise flux received from the planet
        # TODO: to be consistent with get_snr, make it so that bin_width is multiplied elsewhere
        # TODO: add multi-mode kernel option

        if self.data.inst['combination_technique'] == 2:
            flux_planet = np.zeros((transm_eff.shape[0], transm_eff.shape[1]))
            noise_planet = np.zeros((transm_eff.shape[0], transm_eff.shape[1]))
            for m in range(transm_eff.shape[1]):
                flux_planet_iterable = flux_planet_spectrum \
                                       * transm_eff[:, m] \
                                       * integration_time \
                                       * self.data.inst['eff_tot'] \
                                       * self.data.inst['telescope_area'] \
                                       # * self.data.inst['wl_bin_widths'] * 1e6# TODO DO NOT ADD WAVELENGTH BINS
                flux_planet[:, m] = flux_planet_iterable
                noise_planet_iterable = flux_planet_spectrum \
                                        * transm_noise[:, m] \
                                        * integration_time \
                                        * self.data.inst['eff_tot'] \
                                        * self.data.inst['telescope_area'] \
                                        # * self.data.inst['wl_bin_widths'] * 1e6
                noise_planet[:, m] = noise_planet_iterable
        else:
            flux_planet = flux_planet_spectrum \
                          * transm_eff \
                          * integration_time \
                          * self.data.inst['eff_tot'] \
                          * self.data.inst['telescope_area'] \
                          #* self.data.inst['wl_bin_widths']*1e6    # TODO per micrometer, so multiply with 1e6
            noise_planet = flux_planet_spectrum \
                           * transm_noise \
                           * integration_time \
                           * self.data.inst['eff_tot'] \
                           * self.data.inst['telescope_area'] \
                           #* self.data.inst['wl_bin_widths']*1e6

        # calculate the noise from the background sources
        noise_bg_list = self.run_socket(s_name='photon_noise',
                                        method='noise',
                                        index=None)
        if self.data.options.models['instrumental_photon_noise']:
            noise_bg_list.append(self.data.inst['ti_leak'])

        # print('all bg noise', noise_bg_list)
        if type(noise_bg_list) == list:
            noise_bg = np.zeros_like(noise_bg_list[0])
            for _, noise in enumerate(noise_bg_list):
                noise_bg += noise
        else:
            noise_bg = noise_bg_list

        noise_bg = noise_bg * integration_time * self.data.inst['eff_tot']
        # print('noise_bg', noise_bg[22])
        # print('noise_planet', noise_planet[22])
        # Add up the noise and calculate the SNR
        noise = (noise_bg + noise_planet) * 2
        # only once for bracewell!
        if self.data.inst['combination_technique'] == 0:
            noise = noise_bg + noise_planet
        # print('noise as recorded in get spectrum', noise)
        #
        snr_spec = np.sqrt((flux_planet ** 2 / noise)) # * self.data.inst['wl_bin_widths']
        # np.sqrt(
        #     (flux_planet ** 2 / noise).sum())
        print('snr_tot', np.sqrt((flux_planet ** 2 / noise).sum()))
        # TODO how do the returns shape themselves in case of kernel nulling
        # return ([self.data.inst['wl_bins'], snr_spec],
        #         flux_planet,
        #         noise,
        #         noise_bg_list)
        return ([self.data.inst['wl_bins'], snr_spec],
                flux_planet,
                noise,
                noise_bg_list)

        # else:
        #     return ([self.data.inst['wl_bins'], snr_spec],
        #             flux_planet,
        #             [noise, noise_bg_list])

    def preset_configuration(self, case: str):
        """
        Sets preset 2D array configurations. Equal aperture areas. They rely on the baseline and aspect ratio.
        This made it easier to loop all or compare, but does make some of the results a bit misleading. ,
        (one of the weaker aspects of my thesis IMO)
        :param case:
            String value of specific preset config. These are all quite arbitrary

        :return:
        """
        if case == 'emma x':
            x_perturb = 1.0
            y_perturb = 1.0
            nulling_baseline = self.data.inst['bl']  # short baseline (x)
            ratio = self.data.options.array['ratio']  # ratio multiplication for large baseline (y)
            self.data.inst['configuration'] = 0.5 * np.array([[nulling_baseline, nulling_baseline*ratio],
                                                                               [-nulling_baseline, nulling_baseline*ratio],
                                                                               [-nulling_baseline, -nulling_baseline*ratio],
                                                                               [nulling_baseline, -nulling_baseline*ratio]])
        elif case == 'TTN':
            bl = self.data.inst['bl'] * self.data.options.array['ratio']
            rotator = 0
            self.data.inst['configuration'] = np.array([[0.577 * bl * np.cos(rotator), 0.577 * bl * np.sin(rotator)],
                                                        [0.577 * bl * np.cos(2 * np.pi / 3 + rotator),
                                                         0.577 * bl * np.sin(2 * np.pi / 3 + rotator)],
                                                        [0.577 * bl * np.cos(4 * np.pi / 3 + rotator),
                                                         0.577 * bl * np.sin(4 * np.pi / 3 + rotator)]])

        elif case == 'guyon 6':
            self.data.inst['configuration'] = self.data.inst['bl'] * np.array([[0.420, 0.000],
                                                                               [0.486, 0.108],
                                                                               [0.476, 0.487],
                                                                               [-0.585, -0.716],
                                                                               [0.585, 0.716],
                                                                               [0.238, -0.971]])
        elif case == 'pentagon':
            bl = self.data.inst['bl'] * 0.5 * 1.05 * self.data.options.array['ratio']
            rotator = 0
            theta = [0, 2 * np.pi / 5, 4 * np.pi / 5, 6 * np.pi / 5, 8 * np.pi / 5]

            self.data.inst['configuration'] = np.array(
                [[bl * np.cos(rotator), bl * np.sin(rotator)],
                 [bl * np.cos(theta[1] + rotator),
                  bl * np.sin(theta[1] + rotator)],
                 [bl * np.cos(theta[2] + rotator),
                  bl * np.sin(theta[2] + rotator)],
                 [bl * np.cos(theta[3] + rotator),
                  bl * np.sin(theta[3] + rotator)],
                 [bl * np.cos(theta[4] + rotator),
                  bl * np.sin(theta[4] + rotator)]])
            # nulling_baseline = self.data.inst['bl']  # short baseline (x)
            # ratio = self.data.options.array['ratio']  # ratio multiplication for large baseline (y)
            # self.data.inst['configuration'] = 0.5 * np.array([[nulling_baseline, nulling_baseline*ratio],
            #                                                   [0, nulling_baseline*ratio*1.2],
            #                                                   [-nulling_baseline, nulling_baseline*ratio],
            #                                                   [-nulling_baseline*0.8, -nulling_baseline*ratio],
            #                                                   [nulling_baseline*0.8, -nulling_baseline*ratio]])


        elif case == 'angel woolf':
            bl = self.data.inst['bl']
            ratio = self.data.options.array['ratio']
            self.data.inst['configuration'] = np.array([[-0.5 * ratio * bl, 0],
                                                        [-0.5*bl, 0],
                                                        [0.5*bl, 0],
                                                        [0.5*ratio*bl, 0]])

        elif case == 'bracewell':
            x_perturb = 1.0
            y_perturb = 0.0
            self.data.inst['configuration'] = self.data.inst['bl'] * self.data.options.array['ratio'] * np.array([[-0.5 * x_perturb, 0 + y_perturb], [0.5, 0]])

        elif case == 'diamond':
            ratio = self.data.options.array['ratio']
            self.data.inst['configuration'] = self.data.inst['bl'] * np.array([[0.0, 0.5*ratio],
                                                                               [0.5, 0.0],
                                                                               [0.0, -0.5*ratio],
                                                                               [-0.5, 0.0]])
        elif case == 'kite':
            c = self.data.options.array['ratio']
            alpha = 2*np.arctan(1/c)
            theta = [np.pi/2 - alpha, np.pi/2, np.pi/2 + alpha, 3*np.pi/2]
            self.data.inst['configuration'] = np.zeros((len(theta), 2))
            for i, titters in enumerate(theta):
                self.data.inst['configuration'][i, 0] = 0.5 * self.data.inst['bl'] * np.sqrt(1+c**2) * np.cos(titters)
                self.data.inst['configuration'][i, 1] = 0.5 * self.data.inst['bl'] * np.sqrt(1+c**2) * np.sin(titters)

        else:
            print('config preset not recognised, no config set')

    # def random_configuration(self, n_apertures: int, check_mode: bool = False):
    #     """
    #     Creates a random_configuration.
    #
    #     Inputs:
    #
    #     size: int
    #         Defines the amount of inputs.
    #
    #     TBA:
    #     - link to beam combination
    #     - create nulling baseline from combination.
    #
    #     """
    #     rng_radius = np.random.default_rng(self.data.options.other['random_seed'])
    #     rng_angle = np.random.default_rng(self.data.options.other['random_seed']+1)
    #     config_list = []
    #     for i in range(n_apertures):
    #         r = 0.5 * self.data.inst['bl'] * np.sqrt(rng_radius.random())
    #         theta = rng_angle.random() * 2 * np.pi
    #         aperture_coord = [r * np.cos(theta), r * np.sin(theta)]
    #         config_list.append(aperture_coord)
    #
    #     self.data.inst['configuration'] = np.array(config_list)
    #     # print('configuration', np.array(config_list))
    #     if check_mode:
    #         return np.array(config_list)

    def kernel_beam_combiner(self, number_apertures: int = 1, check_for_redundancy: bool = False) -> object:
        """
        This method will open kernel nulling beam combiners. They will follow from research done by Romain Laugier
        His work on arbitrary numbers of apertures will come in handy. How will the performance of the array vary?
        :return:
        """

        # from collections import Counter
        resource_path = Path(__file__).parents[2] / "docs/data"

        matrixfile = resource_path.joinpath(str(number_apertures) + "T_matrices.txt")
        with open(matrixfile, "r") as f:
            combination_matrix = sp.sympify(f.read())
        if check_for_redundancy:
            combination = [tuple(combination_matrix[i]) for i in range(len(combination_matrix))]
            combination_set = set(combination)
            if len(combination_set) < len(combination):
                print('redundant matrices!')
        # fg.plot_outputs_smart(np.array(combination_matrix[i], dtype=np.complex128), onlyoneticklabel=True, title='premade beam combination matrix')

        self.data.inst['U'] = np.array(combination_matrix[0],
                                       dtype=np.complex128)  # np.array(combination_matrix, dtype=np.complex128)

    def beam_combination_setting(self, config: str):
        """
        This method applies preset combination matrices to the instrument settings. the combination matrix is used in
        transmission. Some alternatives are saved as comments in the method.
        :return:
        """

        # Default = shitty recombination matrix
        phi_perturbation = 1.00 # rad
        combination_matrix = (1 / np.sqrt(4)) * np.array([[0, 0, np.sqrt(2), np.sqrt(2)],
                                                          [np.sqrt(2), np.sqrt(2), 0, 0],
                                                          [1, -1, -np.exp(np.pi * 1j / 2 * phi_perturbation), np.exp(np.pi * 1j / 2)],
                                                          [1, -1, np.exp(np.pi * 1j / 2 * phi_perturbation), -np.exp(np.pi * 1j / 2)]])
        # U = (1/2)*np.transpose(np.array([[np.sqrt(2), np.sqrt(2), 0, 0],
        #                              [1, -1, 1j, -1j],
        #                              [1j, -1j, 1, -1],
        #                              [0, 0, np.sqrt(2), np.sqrt(2)]]))

        if config == 'angel woolf':
            # combination_matrix = np.transpose([[0.25, 0.4, 0.25, 0.1],
            #                                    [0.25, 0.1, -0.25, -0.4],
            #                                    [0.25, -0.1, -0.25, 0.4],
            #                                    [0.25, -0.4, 0.25, -0.1]])
            # combination_matrix= np.transpose(np.transpose([[0.2477, 0.3990, 0.2523, 0.1010],
            #                     [0.2523, 0.1010, -0.2477, -0.3990],
            #                     [0.2523, -0.1010, -0.2477, 0.3990],
            #                     [0.2477, -0.3990, 0.2523, -0.1010]]))
            combination_matrix = (1 / 2) * np.array([[0, 0, np.sqrt(2), np.sqrt(2)],
                                                     [np.sqrt(2), np.sqrt(2), 0, 0],
                                                     [1, -1, -1j, 1j],
                                                     [1, -1, 1j, -1j]])

        elif config == 'guyon 6':
            combination_matrix = np.array([[0.1667, 0.0059, 0.0257, 0.2888, 0.0120, 0.5008],
                                               [0.1667, 0.0235, 0.0290, 0.2194, -0.2352, -0.3263],
                                               [0.1667, 0.1232, -0.0183, -0.0014, 0.5800, -0.1104],
                                               [0.1667, -0.3539, -0.4744, 0.0000, -0.0040, -0.0001],
                                               [0.1667, 0.2537, -0.0256, -0.3328, -0.1606, 0.0607],
                                               [0.1667, -0.2399, 0.4270, -0.1575, 0.0072, -0.0018]])

        elif config == 'bracewell':
            combination_matrix = np.array([[1, 1], [-1j, 1j]])

        elif config == 'diamond':
            # combination_matrix = (1 / np.sqrt(4)) * np.array([[1, 0, 1, 0],
            #                                              [0, 1, 0, 1],
            #                                              [-1j, 0, 1j, 0],
            #                                              [0, -1j, 0, 1j]])
            combination_matrix = (1 / np.sqrt(4)) * np.array([[0, 0, np.sqrt(2), np.sqrt(2)],
                                                          [np.sqrt(2), np.sqrt(2), 0, 0],
                                                          [1, -1, -np.exp(np.pi * 1j / 2 * phi_perturbation), np.exp(np.pi * 1j / 2)],
                                                          [1, -1, np.exp(np.pi * 1j / 2 * phi_perturbation), -np.exp(np.pi * 1j / 2)]])

        # elif config == 'pentagon':

            # combination_matrix = (1 / np.sqrt(5)) * np.array([[1, 1, 1, 1, 1],
            #                                                   [1, -np.exp(2* np.pi * 1j / 5), -np.exp(4* np.pi * 1j / 5) , -np.exp(6* np.pi * 1j / 5) , -np.exp(8* np.pi * 1j / 5)   ],
            #                                                   [1,-np.exp(4* np.pi * 1j / 5), -np.exp(8* np.pi * 1j / 5) , -np.exp(2* np.pi * 1j / 5) , -np.exp(6* np.pi * 1j / 5)],
            #                                                   [1,-np.exp(6* np.pi * 1j / 5), -np.exp(2* np.pi * 1j / 5) , -np.exp(8* np.pi * 1j / 5) , -np.exp(4* np.pi * 1j / 5)],
            #                                                   [1, -np.exp(8* np.pi * 1j / 5), -np.exp(6* np.pi * 1j / 5) , -np.exp(4* np.pi * 1j / 5) , -np.exp(2* np.pi * 1j / 5)]])
        return combination_matrix

from typing import Union

import numpy as np

from lifesim.core.modules import PhotonNoiseModule, TransmissionModule
from lifesim.util.radiation import black_body
import lifesim.util.figures as fg


class PhotonNoiseStar(PhotonNoiseModule):
    """
    This class simulates the noise contribution of central star to the interferometric measurement
    of LIFE due to leakage through the null.
    """

    def __init__(self,
                 name: str):
        """
        Parameters
        ----------
        name : str
            Name of the module.
        """

        super().__init__(name=name)
        self.add_socket(s_name='transmission_star',
                        s_type=TransmissionModule)

    def noise(self,
              index: Union[int, type(None)]):
        """
        Simulates the amount of photon noise originating from the star of the observed system
        leaking into the LIFE array measurement.

        Parameters
        ----------
        index: Union[int, type(None)]
            Specifies the planet for which to calculate the noise contribution. If an integer n is
            given, the noise will be calculated for the n-th row in the `data.catalog`. If `None`
            is given, the noise is calculated for the parameters located in `data.single`.

        Returns
        -------
        sl_leak
            Stellar leakage in [photon s-1] per wavelength bin.

        Notes
        -----
        All of the following parameters are needed for the calculation of the exozodi noise
        contribution and should be specified either in `data.catalog` or in `data.single`.

        radius_s : float
            Radius of the observed star in [sun radii].
        distance_s : float
            Distance between the observed star and the LIFE array in [pc].
        temp_s : float
            Temperature of the observed star in [K].
        data.inst['wl_bins'] : np.ndarray
            Central values of the spectral bins in the wavelength regime in [m].
        data.inst['wl_widths'] : np.ndarray
            Widths of the spectral wavelength bins in [m].
        data.inst['telescope_area'] : float
            Area of all array apertures combined in [m^2].

        Raises
        ------
        ValueError
            If the specified transmission map does not exist.
        """

        image_size = 50
        map_selection = 'tm3'

        if index is None:   # for single planet spectrum
            radius_s = self.data.single['radius_s']
            distance_s = self.data.single['distance_s']
            temp_s = self.data.single['temp_s']
        else:
            radius_s = self.data.catalog.radius_s.iloc[index]
            distance_s = self.data.catalog.distance_s.iloc[index]
            temp_s = self.data.catalog.temp_s.iloc[index]

        # check if the specified map exists
        if map_selection not in ['tm1', 'tm2', 'tm3', 'tm4']:
            raise ValueError('Nonexistent transmission map')

        # convert units
        Rs_au = 0.00465047 * radius_s
        Rs_as = Rs_au / distance_s
        Rs_mas = float(Rs_as)  # This is still arcseconds
        Rs_rad = Rs_mas / (3600. * 180.) * np.pi
        # print('size of star in rad', Rs_rad)
        # TODO Instead of recalculating the transmission map for the stellar radius here, one could try
        #   to reuse the inner part of the transmission map already calculated in the get_snr function
        #   of the instrument class
        # TODO: why are we not reusing the maps calculated in the instrument class
        tm_star = self.run_socket(method='transmission_map',
                                  s_name='transmission_star',
                                  hfov=Rs_rad,
                                  image_size=image_size)[1::2]  # odd outputs
        # normalization settings:
        multi = self.data.inst['normalization_multiplier']
        # normalization
        if self.data.inst['combination_technique'] == 2:
            norm_map = 0
        else:
            norm_map = -1
        # tm_star_without = multi * np.array([(tm_star[:, x, :, :] - np.amin(tm_star[:, x, :, :])) /
        #                                     (np.amax(self.data.inst['norm_maps'][norm_map]) - np.amin(
        #                                         tm_star[:, x, :, :]))
        #                                     for x in
        #                                     range(tm_star.shape[1])])
        tm_star = multi * np.array([(tm_star[:, x, :, :] - np.amin(tm_star[:, x, :, :])) /
             (np.amax(self.data.inst['norm_maps'][norm_map]) - np.amin(tm_star[:, x, :, :])) for x in
             range(tm_star.shape[1])]) + self.data.options.models['null_depth']

        # --- null depth pciture
        # map the star
        x_map = np.tile(np.array(range(0, image_size)), (image_size, 1))
        y_map = x_map.T
        r_square_map = (x_map - (image_size - 1) / 2) ** 2 + (y_map - (image_size - 1) / 2) ** 2
        star_px = np.where(r_square_map < (image_size / 2) ** 2, 1, 0)
        bb_star = black_body(bins=self.data.inst['wl_bins'],
                       width=self.data.inst['wl_bin_widths'],
                       temp=temp_s,
                       radius=radius_s,
                       distance=distance_s,
                       mode='star')
        # fig
        # fg.single_map(tm_star[11][-1])

        # print('null depth of stellar map', np.amin(tm_star[11][-1]))
        # print('location of minimum', np.where(tm_star[11][-1] == np.amin(tm_star[11][-1])))
        # print('null depth at centre of stellar map', tm_star[11, -1, 23:28, 23:28])
        # get the stellar leakage
        sl_leak = np.zeros((tm_star.shape[0], tm_star.shape[1]))
        for i, stars in enumerate(np.swapaxes(tm_star, 0, 1)):
            sl_leak[:, i] = (star_px * stars).sum(axis=(-2, -1)) / star_px.sum(
            ) * bb_star * self.data.inst['telescope_area']
        # print('sl_leak', sl_leak)
        if self.data.inst['combination_technique'] == 2:
            pass
        else:
            sl_leak = sl_leak[:, -1]
        if self.data.options.other['talking_star']:
            extent = [-Rs_as, Rs_as, -Rs_as, Rs_as]
            fg.single_map(tm_star[18][-1]*bb_star[11]*star_px, extent=extent,
                          # intensity='ph/m^2/s',
                          plot_title='star map', labels=['alpha [arcsec]', 'beta [arcsec]'], save=True,
                          save_name='Star_map_photon_flux', cmap='Reds', intensity= 'Photon flux [ph/s]')
            print('sl_leak 1', sl_leak[0])
            # ------- figure of
            fg.single_map(tm_star[18][-1], extent=extent, labels=['alpha [arcsec]', 'beta [arcsec]'],
                          plot_title='star null + 10$^{-5}$', cmap='hot', vmin=0, vmax=np.amax(tm_star[18][-1]), save=True, save_name='star map 10_5', intensity='Transmission [-] ')
            fg.single_map(tm_star_without[11][-1], extent=extent, labels=['alpha [arcsec]', 'beta [arcsec]'],
                          cmap='hot', vmin=0, vmax=np.amax(tm_star[18][-1]), save=True, save_name='star map 0', intensity='Transmission [-] ',
                          plot_title='star null')
            import matplotlib.pyplot as plt
            print('shape', np.shape(tm_star[18][-1]))
            xticks = np.arange(-Rs_as, Rs_as, 2 * Rs_as / 50)
            plt.plot(xticks, tm_star[11][-1][24][:], label='null depth $10^{-5}$')
            plt.plot(xticks, tm_star_without[11][-1][24][:], label='null depth 0')
            plt.xlabel('alpha [arcsec]')
            plt.xlim(-Rs_as, Rs_as)
            plt.ylabel('null depth [-]')
            plt.grid()
            plt.legend()
            plt.show()
            print('sl_leak_without', ((star_px * tm_star_without[:, -1]).sum(axis=(-2, -1)) / star_px.sum(
            ) * bb_star * self.data.inst['telescope_area'])[18])
            print('sl_leak with', sl_leak[18])
        return sl_leak


import numpy as np
from typing import Union

from lifesim.core.modules import PhotonNoiseModule
# from lifesim.util import constants
from lifesim.util.figures import single_map
from lifesim.util.radiation import black_body
import matplotlib.pyplot as plt

class PhotonNoiseExozodi(PhotonNoiseModule):
    """
    This class simulates the noise contribution of an exozodi disk to the interferometric
    measurement of LIFE.
    """

    def __init__(self,
                 name: str):
        super().__init__(name=name)
        """
        Parameters
        ----------
        name : str
            Name of the module.
        """

    # def noise(self,
    #            index: Union[int, type(None)]):
    #     pass

    def noise(self,
              index: Union[int, type(None)]):
        """
        Simulates the amount of photon noise originating from the exozodi of the observed system
        leaking into the LIFE array measurement.

        Parameters
        ----------
        index: Union[int, type(None)]
            Specifies the planet for which to calculate the noise contribution. If an integer n is
            given, the noise will be calculated for the n-th row in the `data.catalog`. If `None`
            is given, the noise is calculated for the parameters located in `data.single`.

        Returns
        -------
        ez_leak
            Exozodi leakage in [photon s-1] per wavelength bin.

        Notes
        -----
        All of the following parameters are needed for the calculation of the exozodi noise
        contribution and should be specified either in `data.catalog` or in `data.single`.

        l_sun : float
            Luminosity of the observed star in [solar luminosities].
        distance_s : float
            Distance between the observed star and the LIFE array in [pc].
        z : float
            Zodi level in the observed system in [zodis].
        mas_pix : np.ndarray
            Contains the size of each pixel projected to the sky in [milliarcseconds].
        rad_pix : np.ndarray
            Contains the size of each pixel projected to the sky in [radians].
        data.inst['radius_map'] : np.ndarray
            Contains the distance of a pixel from the center of the detector in [pix].
        data.options.other['image_size']
            Number of pixels on one axis of a square detector (dimensionless). I.e. for a 512x512
            detector this value is 512.
        wl_bins : np.ndarray
            Central values of the spectral bins in the wavelength regime in [m].
        wl_widths : np.ndarray
            Widths of the spectral wavelength bins in [m].
        data.inst['telescope_area'] : float
            Area of all array apertures combined in [m^2].
        data.inst['t_map'] : np.ndarray
            Transmission map of the TM3 mode of the array created by the
            lifesim.TransmissionMap module.
        """

        # read from catalog or single data depending on index specification
        if index is None:
            l_sun = self.data.single['l_sun']
            distance_s = self.data.single['distance_s']
            z = self.data.single['z']
        else:
            l_sun = self.data.catalog.l_sun.iloc[index]
            distance_s = self.data.catalog.distance_s.iloc[index]
            z = self.data.catalog.z.iloc[index]  # How to fix this

        # calculate the parameters required by Kennedy2015
        alpha = 0.34
        r_in = 0.034422617777777775 * np.sqrt(l_sun)
        r_0 = np.sqrt(l_sun)
        sigma_zero = 7.11889e-8  # Sigma_{m,0} from Kennedy+2015 (doi:10.1088/0067-0049/216/2/23)

        # reshape the mas per pixel array for calculation (to (n, 1, 1))
        mas_pix = np.array([self.data.inst['mas_pix']])
        if mas_pix.shape[-1] > 1:
            mas_pix = np.reshape(mas_pix, (mas_pix.shape[-1], 1, 1))
        rad_pix = np.array([self.data.inst['rad_pix']])
        if rad_pix.shape[-1] > 1:
            rad_pix = np.reshape(rad_pix, (rad_pix.shape[-1], 1, 1))

        au_pix = mas_pix / 1e3 * distance_s

        # the radius as measured from the central star for every pixel in [AU]
        r_au = self.data.inst['radius_map'] * au_pix

        # identify all pixels where the radius is larges than the inner radius by Kennedy+2015
        r_cond = ((r_au >= r_in)
                  & (r_au <= self.data.options.other['image_size'] / 2 * au_pix))

        # calculate the temperature at all pixel positions according to Kennedy2015 Eq. 2
        temp_map = np.where(r_cond,
                            278.3 * (l_sun ** 0.25) / np.sqrt(r_au), 0)
        # Temp map figure
        # print('what are we dealing with here', l_sun, r_in, r_0)
        # print('exo temp map shape', np.shape(temp_map))
        # print('r_au', np.shape(r_au))
        # print('crosssections temp', r_au[18, :, 128], temp_map[18, :, 128])
        # fig, ax = plt.subplots()
        # ax.loglog(r_au[18, :, 128], temp_map[18, :, 128])
        # ax.set_xlabel('Separation from star [AU]', fontsize=18)
        # ax.set_yticks([100, 1000])
        # ax.set_xticks([0.1, 1.0])
        # ax.set_xlim([0.02, 9])
        # ax.set_ylim([200, 1300])
        # ax.set_ylabel('Temperature [K]', fontsize=18)
        # fig.tight_layout()
        # plt.show()


        # print('r_cond', sum(r_cond), 'sigma_zero', sigma_zero, 'zodi', z, 'r_AU', r_au, 'r_0', r_0, ' alpha', alpha)
        # calculate the Sigma (Eq. 3) in Kennedy2015 and set everything inside the inner radius to 0
        sigma = np.where(r_cond,
                         sigma_zero * z *
                         (r_au / r_0) ** (-alpha), 0)
        # sigma figure
        # fig, ax = plt.subplots()
        # ax.loglog(r_au[18, :, 128], sigma[18, :, 128])
        # print('crosssections sig', r_au[18, :, 128], sigma[18, :, 128])
        # ax.set_xlabel('Separation from star [AU]', fontsize=18)
        # # ax.set_yticks([100, 1000])
        # ax.set_xticks([0.1, 1.0])
        # ax.set_xlim([0.04, 9])
        # ax.set_ylim([6e-8, 2e-7])
        # ax.set_ylabel('Surface denisty [AU$^2$/AU$^2$]', fontsize=18)
        # fig.tight_layout()
        # plt.show()

        wl_bins = np.array([self.data.inst['wl_bins']])
        if wl_bins.shape[-1] > 1:
            wl_bins = np.reshape(wl_bins, (wl_bins.shape[-1], 1, 1))

        wl_bin_widths = np.array([self.data.inst['wl_bin_widths']])
        if wl_bin_widths.shape[-1] > 1:
            wl_bin_widths = np.reshape(wl_bin_widths, (wl_bin_widths.shape[-1], 1, 1))
        # print('noise area (exozodi)', self.data.inst['telescope_area'])
        # print('radpix (exozodi)', self.data.inst['rad_pix'])
        # get the black body radiation emitted by the interexoplanetary dust
        # print('rad_pix', rad_pix)
        f_nu_disk = black_body(bins=wl_bins,
                               width=wl_bin_widths,
                               temp=temp_map,
                               mode='wavelength') \
                    * sigma * rad_pix ** 2 * self.data.inst['telescope_area']

        # print('disk shape', np.shape(f_nu_disk))

        # brightness figure
        # fig, ax = plt.subplots()
        # print('radpix', rad_pix ** 2)
        # print('radpix', mas_pix ** 2)
        # ax.loglog(r_au[18, :, 128], f_nu_disk[18, :, 128])
        # print('crosssections flux', r_au[18, :, 128], sigma[18, :, 128])
        # ax.set_xlabel('Separation from star [AU]', fontsize=18)
        # # ax.set_yticks([100, 1000])
        # ax.set_xticks([0.1, 1.0])
        # ax.set_xlim([0.04, 2.5])
        # # ax.set_ylim([6e-8, 2e-7])
        # ax.set_ylabel('Exozodi flux [ph/s/sr]', fontsize=18)
        # fig.tight_layout()
        # plt.show()

        # print('f_nu', np.shape(f_nu_disk), 'fnu poart', sum(f_nu_disk[0]))
        ap = np.where(self.data.inst['radius_map']
                      <= self.data.options.other['image_size'] / 2, 1, 0)
        # figures
        # wl_slab = 18
        # extent_value = self.data.inst['hfov_mas'][wl_slab]
        # extent = [-extent_value, extent_value, -extent_value, extent_value]
        # labels = ['alpha [mas]', 'beta[mas]']
        # single_map(f_nu_disk[wl_slab], extent=extent, plot_title=None, labels=labels,
        #            save=True, save_name='Exodisc flux',
        #            cmap='hot', intensity='ph/s [ph/s]')
        # single_map(self.data.inst['maps'][11, -1] * ap, extent=extent, labels=labels, save=True,
        #            save_name='exozodi emma map', cmap='hot', intensity='Normalized transmission[-] ')
        # single_map(f_nu_disk[wl_slab] * ap * (self.data.inst['maps'][11, -1]), extent=extent, plot_title=format(self.data.inst['wl_bins'][wl_slab], '.2f')+ '$\mu$m', labels=labels, save=True,
        #            save_name='Emma & exozodi convolution sr', cmap='hot',
        #            intensity='Spectral radiance [ph/s] ')


        ez_leak = np.zeros((self.data.inst['maps'].shape[0], self.data.inst['maps'][:, 1::2].shape[1]))
        for i, iter_map in enumerate(np.swapaxes(self.data.inst['maps'][:, 1::2, :, :], 0, 1)):
            # single_map(iter_map[11], plot_title='transmission exozodiacal noise map '+str(i))
            # print(iter_map.shape)
            # single_map(iter_map[11])
            # single_map(f_nu_disk[11])
            ez_leak[:, i] = (iter_map * f_nu_disk * ap).sum(axis=(-2, -1))
        if self.data.inst['combination_technique'] == 2:
            pass
        else:
            ez_leak = ez_leak[:, -1]
            # print('ez_leak1', ez_leak)
            # ez_leak = (f_nu_disk * self.data.inst['maps'][:, -2] * ap).sum(axis=(-2, -1))
        # print('ez_leak_shape', ez_leak.shape)
        # print('ez_leak', ez_leak)
        return ez_leak

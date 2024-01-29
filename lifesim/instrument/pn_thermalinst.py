
import numpy as np

from lifesim.core.modules import PhotonInstrumentNoiseModule
from lifesim.util.radiation import black_body
from lifesim.util.figures import single_map


class PhotonNoiseThermalInst(PhotonInstrumentNoiseModule):
    """
    This class simulates the noise contribution of the thermal localzodical dust to the
    interferometric measurement of LIFE.
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

    def inst_t_noise(self):
        """
        Simulates the amount of photon noise originating from the instrument leaking into the LIFE
        optical train. The instrument is simulated as a grey body with emissivity 0.25 à la colin.

        Parameters
        ----------
        index: Union[int, type(None)]
            Specifies the planet for which to calculate the noise contribution. If an integer n is
            given, the noise will be calculated for the n-th row in the `data.catalog`. If `None`
            is given, the noise is calculated for the parameters located in `data.single`.

        Returns
        -------
        ti_leak
            thermal instrument leakage in [photon s-1] per wavelength bin.

        Notes
        -----
        All of the following parameters are needed for the calculation of the thermal instrument noise
        contribution and should be specified either in `data.catalog` or in `data.single`.

        data.options.array['instrument_temperature'] : float
            Designates the instrument temperature
        data.inst['radius_map'] : np.ndarray
            Contains the distance of a pixel from the center of the detector in [pix].
        data.options.other['image_size']
            Number of pixels on one axis of a square detector (dimensionless). I.e. for a 512x512
            detector this value is 512.
        data.inst['wl_bins'] : np.ndarray
            Central values of the spectral bins in the wavelength regime in [m].
        data.inst['wl_widths'] : np.ndarray
            Widths of the spectral wavelength bins in [m].
        data.inst['hfov'] : np.ndarray
            Contains the half field of view of the observatory in [rad] for each of the spectral
            bins.
        data.inst['t_map'] : np.ndarray
            Transmission map of the TM3 mode of the array created by the
            lifesim.TransmissionMap module.
        data.inst['telescope_area'] : float
            Area of all array apertures combined in [m^2].

        Raises
        ------
        ValueError
            If the specified localzodi model does not exist.
        """

        ap = np.where(self.data.inst['radius_map']
                      <= self.data.options.other['image_size'] / 2, 1, 0)

        epsilon = 0.25
        # epsilon = self.data.options.array['instrument_emissivity']
        ti_flux_sr = epsilon * black_body(mode='wavelength',
                                          bins=self.data.inst['wl_bins'],
                                          width=self.data.inst['wl_bin_widths'],
                                          temp=self.data.options.array['instrument_temperature'])
        # from steridian to absolute per area ?
        ti_flux = ti_flux_sr * (np.pi * self.data.inst['hfov'] ** 2)
        # ----------- figures
        # import matplotlib.pyplot as plt
        # fig, ax1 = plt.subplots()
        # # ax1.stairs(ti_flux_sr, self.data.inst['wl_bin_edges'] * 10 ** 6, color='k', label='Spectral radiance')
        # ax1.stairs(ti_flux* self.data.inst['telescope_area'], self.data.inst['wl_bin_edges'] * 10 ** 6, color='r', label='Flux_s before mask')
        # ax1.set_yscale('log')
        # ax1.set_ylim([10**-1, 10**10])
        # ax1.set_ylabel('Spectral signal [ph/s/m]')
        # ax1.legend(loc='upper right')
        # fig
        # print('ap divided by ap sum', ap/ap.sum())
        # # single_map(self.data.inst['maps'][11, -1])
        # # single_map(ti_flux[27]*ap*self.data.inst['maps'][27, -1], cmap='Reds')
        # wl_slab=int(len(self.data.inst['hfov_mas'])*0.75)
        # print('slab', self.data.inst['wl_bins'][wl_slab])
        # # print('slab', wl_slab)
        # extent_value = self.data.inst['hfov_mas'][wl_slab]
        # # print(np.shape(self.data.inst['hfov_mas']))
        # extent = [-extent_value, extent_value, -extent_value, extent_value]
        # labels = ['alpha [mas]', 'beta[mas]']
        # print('ti_flux before map', ti_flux[wl_slab])
        # # single_map(ti_flux[wl_slab] * ap /(256*256), extent=extent, labels=labels, intensity='Spectral flux density [ph/s/m$^2$/m] ', cmap='Reds', save=True, save_name='cropped instrument noise flux')
        # print('ti_flux after map', (ti_flux[wl_slab] * ap).sum(axis=(-2,-1)))
        # -------------------- Run
        ti_leak = np.zeros((self.data.inst['maps'].shape[0], self.data.inst['maps'][:, 1::2].shape[1]))
        for i, iter_map in enumerate(np.swapaxes(self.data.inst['maps'][:, 1::2, :, :], 0, 1)):
            # print('flux before bro', ti_flux* self.data.inst['telescope_area'])
            # ti_leak[:, i] = ti_flux*0.785398163397\
            #       * self.data.inst['telescope_area']

            ti_leak[:, i] = (ap * iter_map).sum(axis=(-2, -1)) / ap.sum() * ti_flux \
                      * self.data.inst['telescope_area']
            # print('what is it', (ap * iter_map).sum(axis=(-2, -1)) / ap.sum()  )
            # print(ti_leak[:, i])

        if self.data.inst['combination_technique'] == 2:
            pass
        else:
            ti_leak = ti_leak[:, -1]

        return ti_leak

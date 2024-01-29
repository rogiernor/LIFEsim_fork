import numpy as np
from typing import Union

from lifesim.core.modules import TransmissionModule
import lifesim.util.figures as fg
from lifesim.util.normalizer import transmission_normalizer, onedimensional_map_normalizer
import matplotlib.pyplot as plt
import scipy as sp

class TransmissionMap(TransmissionModule):
    """
    Module for calculating nulling interferometer transmission maps.
    """
    def __init__(self,
                 name: str):
        """
        Parameters
        ----------
        name : str
            Name of the module.
        """

        super().__init__(name)

    def complex_amplitude_vector(self,
                                 amplitude,
                                 wavelength,
                                 alpha,
                                 beta,
                                 x,
                                 y):
        """
        Method for calculating the complex amplitude vector, originally after guyon et al. 2013
        V is calculated for each aperture for each wavelength, see output.
        Parameters:
        amplitude: any
            The amplitude of the beam when normalized, or the radius of the constituent telescope for varying radii,
            but that remains untested. # TODO For future research
        wavelength: array
            wavelengths of bins
        alpha: array
            latitude in FOV
        beta: array
            longitude in FOV
        x: float
            x-location of given constituent telescope
        y: float
            y-location of given consituent telescope

        :return:
        Complex Amplitude Vector V_k, which is a 3 dimension ndarray: [lambda, alpha, beta]
        """

        vk = amplitude * \
            np.exp(1j*2*np.pi*(x * alpha + y * beta)
                   / wavelength)
        # fg.single_map(abs(vk)[:,:,18], plot_title='for x-loc'+str(x)+'and y-loc'+str(y))
        return vk

    def transmission_map(self,
                          direct_mode: bool = False,
                          d_alpha: np.ndarray = None,
                          d_beta: np.ndarray = None,
                          hfov: np.ndarray = None,
                          image_size: int = None):
        """
        Return the transmission intensity outputs of a designated configuration for a nulling
        interferometry telescope array. Made by Rogier

        Parameters
        ----------
        direct_mode: bool
            If direct mode is set to True, the variables wl_bins, d_alpha and d_beta are injected
            into the transmission calculation directly. The output is not necessarily in the
            form of a map. The inputs hfov and image size will be ignored in this mode.
        d_alpha: np.ndarray
            The x-positions of the points to be evaluated measured from the central viewing axis in
            [rad].
        d_beta: np.ndarray
            The y-positions of the points to be evaluated measured from the central viewing axis in
            [rad].
        hfov : np.ndarray
            Contains the half field of view of the observatory in [rad] for each of the spectral
            bins. If no value is given, `data.inst['hfov']` is used.
        image_size : int
            Number of pixels on one axis of a square detector (dimensionless). I.e. for a 512x512
            detector this value is 512. If no value is given, `data.options.other['image_size']` is
            used.

        Returns
        -------
        tm
            4-dimensional array containing all modes' transmission maps


        Notes
        -----
        The following additional parameters are required for the calculation of the transmission
        maps and should be specified either in `data.catalog` or in `data.single`.

        data.inst['wl_bins'] : np.ndarray
            Central values of the spectral bins in the wavelength regime in [m].
        data.inst['configuration'] : np.ndarray
            X and y coordinates of the array apertures.
        data.inst['diameters']
            EXPERIMENTAL size of apertures, usually assigned equal. Can be further edited in case of more true to life
            amplitude considerations.
        data.options.array['n_apertures'] : int
            Number of apertures. This is both ways dependent on the type of beam combiner used.
        self.data.inst['U'] : np.ndarray
            Beam combination, in matrix form.

        """

        if hfov is None:
            hfov = self.data.inst['hfov']
        if image_size is None:
            image_size = self.data.options.other['image_size']

        # reshape the wl_bins and hfov arrays for calculation (to (n, 1, 1))
        wl_bins = np.array([self.data.inst['wl_bins']])  # wavelength in m
        if wl_bins.shape[-1] > 1:
            wl_bins = np.reshape(wl_bins, (wl_bins.shape[-1], 1, 1))

        if direct_mode:
            alpha = d_alpha
            beta = d_beta
            # print('length of alpha', alpha, len(alpha), 'and beta ', len(beta), 'in direct_mode')
            vector_size = [1, len(alpha)]
            # print('vector size', vector_size)
        else:
            hfov = np.array([hfov])  # wavelength in m
            if hfov.shape[-1] > 1:
                hfov = np.reshape(hfov, (hfov.shape[-1], 1, 1))

            # generate 1D array that spans field of view
            angle = np.linspace(-1, 1, image_size)

            # angle matrix in x-direction ("alpha")
            alpha = np.tile(angle, (image_size, 1))

            # angle matrix in y-direction ("beta")
            beta = alpha.T

            # convert angle matrices to fov units
            alpha = alpha * hfov
            beta = beta * hfov
            #sizes
            vector_size = [len(alpha[0]), len(beta[0])]

        # config = 'emma x'
        U = self.data.inst['U']
        # initiate loop
        tm = np.empty([len(U[:, 1]), len(wl_bins), vector_size[0], vector_size[1]])
        for n in range(U.shape[0]):  # loop over rows
            aperture_beam = []
            for m in range(self.data.options.array['n_apertures']):  # loop over columns  (aperture) TODO add better dimensioning for apertures
                # print('at column (aperture)', m,
                #       ', x is:', self.data.inst['configuration'][m][0],
                #       ', y is:', self.data.inst['configuration'][m][1])
                v = self.complex_amplitude_vector(amplitude=self.data.inst['diameters'][m]/2,
                                                  wavelength=wl_bins,
                                                  alpha=alpha,
                                                  beta=beta,
                                                  x=self.data.inst['configuration'][m][0],
                                                  y=self.data.inst['configuration'][m][1])
                coupling_efficiency = 1.0
                aperture_beam.append(U[n, m] * v * coupling_efficiency)
            tm[n, :, :, :] = np.square(abs(sum(aperture_beam)))
        return tm

    def transmission_efficiency(self,
                                index: Union[int, type(None)]):
        """
        Integrates over transmission curves to get the transmission efficiency for signal and
        noise.

        Parameters
        ----------
        index: Union[int, type(None)]
            Specifies the planet for which to calculate the transmission efficiency. If an integer
            n is given, the noise will be calculated for the n-th row in the `data.catalog`. If
            `None` is given, the noise is calculated for the parameters located in `data.single`.

        Returns
        -------
        transm_eff
            Transmission efficiency per spectral bin for the exoplanet signal
        transm_noise
            Transmission efficiency per spectral bin for the photon noise received from the
            exoplanet signal

        Notes
        -----
        The following additional parameters are required for the calculation of the transmission
        efficiency and should be specified either in `data.catalog` or in `data.single`.

        data.single['angsep']
            Angular separation between the observed star and the observed exoplanet in [arcsec].
        data.inst['wl_bins'] : np.ndarray
            Central values of the spectral bins in the wavelength regime in [m].
        """

        if index is None:
            angsep = self.data.single['angsep']
        else:
            angsep = self.data.catalog.angsep.iloc[index]

        # Rogier version
        tc_chop, tc_tm4 = self.transmission_curve(angsep=angsep)
        # integrate over angles to get transmission efficiency
        transm_eff = np.sqrt((tc_chop ** 2).mean(axis=(-2, -1)))
        transm_noise = np.sqrt((tc_tm4 ** 2).mean(axis=(-2, -1)))
        return transm_eff, transm_noise

    def transmission_curve(self,
                            angsep: float,
                            phi_n: int = 360,
                            n_cycles: float = 2*np.pi):
        """
        Calculates the radial transmission curve of the LIFE array. It is therefore key to the performance and speed
        of the get_snr() operation.

        Parameters
        ----------
        angsep : float
            Angular separation between the observed star and the observed exoplanet in [arcsec].
        phi_n : int
            Number of rotation steps used in integration.
        use_in_snr: bool = True
            Extra option. If the curve is called before the creation of an un-normalized map, this function can be set
            to False, and the curve will generate the map itself.

        Returns
        -------
        transm_curve_chop
            Radial transmission curve corresponding to the chopped transmission map.
        transm_curve_tm4
            Radial transmission curve corresponding to the transmission map of the 4th mode 'tm4'.
        """

        # convert angular separation to radians
        angsep_rad = angsep / (3600 * 180) * np.pi

        # create 1D array with azimuthal coordinates
        phi_lin = np.linspace(0, n_cycles, phi_n, endpoint=False)

        # retrieve the transmission curves
        transmission_curves = self.transmission_map(direct_mode=True,
                                                     d_alpha=angsep_rad * np.cos(phi_lin),
                                                     d_beta=angsep_rad * np.sin(phi_lin))
        # import matplotlib.pyplot as plt
        # normalizer settings:
        multi = self.data.inst['normalization_multiplier']
        technique = self.data.inst['combination_technique']
        if technique == 2:
            # print('curves in curves', transmission_curves)
            # plt.plot(phi_lin, transmission_curves[1][11][0])
            # plt.title('mode 1 at ' + str(angsep) + ' arcsec')
            # plt.show()
            transmission_curves = multi * transmission_normalizer(transmission_curves,
                                                                  normalization_map=self.data.inst['norm_maps'][0])
            # plt.plot(phi_lin, transmission_curves[11][1][0])
            # plt.title('mode 1 at ' + str(angsep) + ' arcsec')
            # plt.show()
            transmission_curves_noise = transmission_curves[:, 1::2, :, :]
            # plt.plot(phi_lin, transmission_curves_noise[18, 0, 0], label='eff')
            # plt.plot(phi_lin, transmission_curves[18, 2, 0, :], label='noise')
            # plt.title('noise eff mode 1 at ' + str(angsep) + ' arcsec')
            # plt.legend()
            # plt.show()
            transmission_curves_signal = [transmission_curves[:, i + 1] - transmission_curves[:, i + 2]
                                          for i in range(0, transmission_curves.shape[1] - 2, 2)]
            transmission_curves_signal = np.swapaxes(np.array(transmission_curves_signal), 0, 1) #TODO this is the only place where normalization maps are switched by axes, and only because theyre switched back
        elif technique == 1:
            transmission_curves = multi * transmission_normalizer(transmission_curves,
                                                                  normalization_map=self.data.inst['norm_maps'][0]*2)

            transmission_curves_signal = transmission_curves[:, 2] - transmission_curves[:, 3]
            transmission_curves_noise = transmission_curves[:, 2, :, :]
            # plt.plot(phi_lin, transmission_curves_signal[0][0])
            # plt.title('X_array 1 at ' + str(angsep) + ' arcsec')
            # plt.show()
        else:
            transmission_curves = multi * transmission_normalizer(transmission_curves,
                                                                  normalization_map=self.data.inst['norm_maps'][0])
            # plt.plot(phi_lin, transmission_curves[11][-1][0])
            # plt.title('X_array 1 at ' + str(angsep) + ' arcsec')
            # plt.show()
            transmission_curves_signal = transmission_curves[:, 1]
            transmission_curves_noise = transmission_curves_signal

        # print('transmission curves shape', transmission_curves.shape)
        return transmission_curves_signal, transmission_curves_noise

    def plot_psf(self, slab: int, save: bool = False, save_name: str = 'RMS point spread function',
                 sliced: bool = False, custom: bool = False):
        '''
        This function is used to create a figure of the RMS modulation efficiency, after Lay + TPF-I (2007).
        :param slab:
        The wavelength bin for which the PSF is shown (required).
        :param save:
        A bool defining whether the figure will be saved (optional: set to False).
        :param save_name:
        Save name under which the plot will appear in 'results' (optional: set to RMS point spread function).
        :param sliced:
        Bool defining whether a sideways cut of the pie chart is also shown. This can be a more useful curve in some
        cases.
        :param custom:
        If true, a custom planet can be evaluated with the plot-psf, to show the modulation efficiency for that planet
        for given config & map. 
        :return:
        Plot of PSF.
        '''
        # ---- extra imports
        # import matplotlib.pyplot as plt
        # ---- initiate axes
        angsep = np.linspace(0, self.data.inst['hfov_mas'][slab], 128, endpoint=False)
        rms_es = np.zeros(len(angsep))
        # ---- initiate figure
        fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
        # ----- loop RMS
        for i, separation in enumerate(angsep):
            sep_as = separation * 0.001     # convert to arcsecond
            curve, _ = self.transmission_curve(sep_as)
            curve_rms = np.sqrt((curve[slab, 0] ** 2).mean())
            rms_es[i] = curve_rms
        # ------ define
        azm = np.linspace(0, 2 * np.pi)
        r, th = np.meshgrid(angsep, azm)
        z = np.tile(rms_es, (r.shape[0], 1))
        ax.set_xlabel('Radial: distance from optical axis [mas]')
        print('r', np.shape(r))
        print('th', np.shape(th))
        print('z', np.shape(z))
        plt.pcolormesh(th, r, z, shading='auto')
        plt.colorbar(label='RMS [-]')
        if save:
            plt.savefig('C:/Users/rogie/Desktop/Afstudeerproject/05. Figures (made by me)/02. Thesis/'
                        + save_name)
        plt.show()
        if sliced:
            plt.plot(angsep, rms_es)
            plt.ylabel('modulation $\eta$ rms [-]', fontsize=18)
            plt.xlabel('angsep [mas]', fontsize=18)
            # plt.title('slice of the pie')
            if save:
                plt.savefig('C:/Users/rogie/Desktop/Afstudeerproject/05. Figures (made by me)/02. Thesis/'
                            + save_name + ' slice of pie')
            plt.show()
            if custom:
                fig, (ax1, ax2) = plt.subplots(1, 2)
                # ----- Transmission curve for the given SINGLE planet in catalog (doesnt work otherwise)
                print('Angsep ', self.data.catalog.angsep.values)
                as_sep = self.data.catalog.angsep.values
                curvo, evil_curvo = self.transmission_curve(as_sep)
                print('rms of curve at angsep', np.sqrt((curvo ** 2).mean(axis=(-2, -1))))
                ramses = np.sqrt((curvo ** 2).mean(axis=(-2, -1)))
                if self.data.options.array['technique'] ==1:
                    ax1.plot(np.linspace(0, 2 * np.pi, 360, endpoint=False), curvo[slab, 0], label='Eff')
                    ax1.plot(np.linspace(0, 2 * np.pi, 360, endpoint=False), evil_curvo[slab, 0], label='noise Eff')
                    ax2.plot(self.data.inst['wl_bins'], ramses, label='rms modulation eff')
                if self.data.options.array['technique'] == 2:
                    for kernel in range(len(curvo[slab,0, :])):
                        # print('shape of curvo', np.shape(curvo))
                        ax1.plot(np.linspace(0, 2 * np.pi, 360, endpoint=False), curvo[slab, kernel, 0], label='Eff output k='+str(kernel))
                        ax1.plot(np.linspace(0, 2 * np.pi, 360, endpoint=False), evil_curvo[slab, kernel, 0],
                                 label='noise eff k=' + str(kernel))
                        # print('shape of ramses', np.shape(ramses))
                        ax2.plot(self.data.inst['wl_bins'], ramses[:, kernel], label='rms modulation eff'+str(kernel))
                ax1.set_ylabel('modulation [-]')
                ax1.set_label('$\phi$ [rad]')
                ax1.legend()
                ax2.legend()
                ax1.grid()
                if save:
                    plt.savefig('C:/Users/rogie/Desktop/Afstudeerproject/05. Figures (made by me)/02. Thesis/'
                                + save_name + ' rotation at'+ str(int(self.data.catalog.angsep*1000))+' mas.png')
                plt.show()
            else:
                index = np.where(rms_es == max(rms_es))[0][0]
                print('max of rms', max(rms_es))
                print('index', index)
                print('angsep = ', angsep[index])
                max_curve, _ = self.transmission_curve(angsep[index]*0.001)
                if self.data.options.array['technique'] ==1:
                    plt.plot(np.linspace(0, 2 * np.pi, 360, endpoint=False), max_curve[slab, 0])
                if self.data.options.array['technique'] == 2:
                    for kernel in range(len(max_curve[slab, 0, :])):
                        plt.plot(np.linspace(0, 2 * np.pi, 360, endpoint=False), max_curve[slab, kernel, 0],
                                 label='Eff output k=' + str(kernel))
                plt.ylabel('modulation [-]')
                plt.xlabel('$\phi$ [rad]')
                plt.grid()
                # plt.title('slice of the pie')
                if save:
                    plt.savefig('C:/Users/rogie/Desktop/Afstudeerproject/05. Figures (made by me)/02. Thesis/'
                                + save_name + ' max rotation.png')
                plt.show()


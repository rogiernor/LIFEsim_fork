
# TODO Pre-kernel transmission curve (verified against Felix)
#
# def transmission_curve2(self,
#                         angsep: float,
#                         phi_n: int = 360,
#                         use_in_snr: bool = True):
#     """
#     Calculates the radial transmission curve of the LIFE array
#
#     Parameters
#     ----------
#     angsep : float
#         Angular separation between the observed star and the observed exoplanet in [arcsec].
#     phi_n : int
#         Number of rotation steps used in integration.
#     use_in_snr: bool = True
#         Extra option. If the curve is called before the creation of an un-normalized map, this function can be set
#         to False, and the curve will generate ther map itself.
#
#     Returns
#     -------
#     transm_curve_chop
#         Radial transmission curve corresponding to the chopped transmission map.
#     transm_curve_tm4
#         Radial transmission curve corresponding to the transmission map of the 4th mode 'tm4'.
#     """
#
#     # convert angular separation to radians
#     angsep_rad = angsep / (3600 * 180) * np.pi
#
#     # create 1D array with azimuthal coordinates
#     phi_lin = np.linspace(0, 2 * np.pi, phi_n, endpoint=False)
#
#     # retrieve the transmission curves
#     transmission_curves = self.transmission_map2(config=self.data.options.other['scenario'],
#                                                  direct_mode=True,
#                                                  d_alpha=angsep_rad * np.cos(phi_lin),
#                                                  d_beta=angsep_rad * np.sin(phi_lin))
#     # check if the method is used in get_snr
#     chop = transmission_curves[-2] - transmission_curves[-1]
#     if use_in_snr:
#         # normalizer settings:
#         multi = self.data.inst['normalization_multiplier']
#         # Planetary output for noise
#         transm_curve_tm4 = multi * self.transmission_normalizer(transmission_curves[-1],
#                                                                 normalization_map=self.data.inst['norm_maps'][[0]])
#         transm_curve_tm3 = multi * self.transmission_normalizer(transmission_curves[-2],
#                                                                 normalization_map=self.data.inst['norm_maps'][[0]])
#         # Planetary output for signal (chopped) - > if more than one chop
#         transm_curve_chop = transm_curve_tm3 - transm_curve_tm4
#     else:
#         # manually normalize
#         normalization_curves = self.transmission_map2(config=self.data.options.other['scenario'])
#         normalization_chop = normalization_curves[-2] - normalization_curves[-1]
#         transm_curve_tm4 = np.array([(transmission_curves[-1, x, :, :] - np.amin(normalization_curves[-1, x, :, :])) / (
#                 np.amax(normalization_curves[-1, x, :, :]) - np.amin(normalization_curves[-1, x, :, :])) for x in
#                                      range(len(transmission_curves[-1]))])
#         transm_curve_chop = np.array([(2 * (chop[x, :, :] - np.amin(normalization_chop[x, :, :]))) / (
#                 np.amax(normalization_chop[x, :, :]) - np.amin(normalization_chop[x, :, :])) - 1 for x in
#                                       range(len(chop))])
#     return transm_curve_chop, transm_curve_tm4\


# ------------------------------- Kernel generator

# def kernel_generator(self, n_a: int):
#     '''
#     Generates kernel nulling arrays. Was used before Romain told me to better just use his matrices.
#     '''
#     # Nbl = np.math.factorial(Na) / \
#     #       (np.math.factorial(2) * np.math.factorial(Na - 2))
#     # Nthk = Nbl - (Na - 1)
#     # thetas = sp.symbols("theta")
#     # The phase shifts according to Na
#     thethetas = [k * (2 * np.pi / n_a) for k in range(n_a)]
#     print('thethetas', thethetas)
#     # normalization multiplier
#     amp = 1 / np.sqrt(n_a)
#     # Possible combinatory solutions
#     perms = np.array(list(itertools.permutations(np.arange(1, n_a, 1, dtype=np.int32))))
#     perms = np.hstack((np.zeros((perms.shape[0], 1), dtype=np.int64), perms))
#     # print(perms)
#     # Adding zero phasors for the bright output in row 1
#     perms = np.vstack((np.zeros((1, perms.shape[1]), dtype=np.int64), perms))
#     phasors = []
#     for i in range(perms.shape[0]):
#         phasors.append([amp * np.exp(1j * thethetas[perms[i, j]]) for j in range(perms.shape[1])])
#     phasors = np.array(phasors)
#     self.other['U'] = phasors
#     print(' shape of combination matrix', self.other['U'].shape)
#     print('generated beam combination matrix', self.other['U'])
#     fg.plot_outputs_smart(self.other['U'], onlyoneticklabel=True, title='generated beam combination matrix')


# ---------------------------------- Adjust baseline to HZ:

# def adjust_bl_to_hz(self,
    #                     hz_center: float,
    #                     distance_s: float):
    #     """
    #     Adjusts the baseline of the array to be optimal for observations in the habitable zone of
    #     the target star for the selected optimal wavelength.
    #
    #     Parameters
    #     ----------
    #     hz_center : float
    #         Separation of the center of the habitable zone in [AU].
    #     distance_s : float
    #         Distance between the observed star and the LIFE array in [pc].
    #     """
    #
    #     # convert the habitable zone to radians
    #     hz_center_rad = hz_center / distance_s / (3600 * 180) * np.pi  # in rad
    #
    #     # put first transmission peak of optimal wl on center of HZ
    #     # for the origin of the value 0.5.. see Ottiger+2021
    #     baseline = (0.589645 / hz_center_rad
    #                             * self.data.options.other['wl_optimal'] * 10 ** (-6))
    #
    #     self.apply_baseline(baseline=baseline)

# ---------------------------------- Apply baseline

    # def apply_baseline(self,
    #                    baseline: float,
    #                    print_warning: bool = False):
    #     """
    #     Adjusts the nulling baseline of the array to the specified value.
    #
    #     Parameters
    #     ----------
    #     baseline : float
    #         Length of the nulling baseline in [m].
    #     print_warning : bool
    #         If set to true, function will print a warning if the specified baseline lies outside
    #         the allow baseline range.
    #     """
    #     # make sure that the baseline does not exeed the set baseline limits
    #     self.data.inst['bl'] = np.maximum(baseline,
    #                                       self.data.options.array['bl_min'])
    #     self.data.inst['bl'] = np.minimum(baseline,
    #                                       self.data.options.array['bl_max'])
    #     if (self.data.inst['bl'] != baseline) and print_warning:
    #         warn('Specified baseline exceeded baseline limits. Baseline fixed to '
    #              'respective limit')
    #
    #     # update the position of the apertures
    #     # self.data.inst['apertures'] = np.array([
    #     #     [-self.data.inst['bl'] / 2,
    #     #      -self.data.options.array['ratio'] * self.data.inst['bl'] / 2., 1.],
    #     #     [self.data.inst['bl'] / 2,
    #     #      -self.data.options.array['ratio'] * self.data.inst['bl'] / 2., 1.],
    #     #     [self.data.inst['bl'] / 2,
    #     #      self.data.options.array['ratio'] * self.data.inst['bl'] / 2., 1.],
    #     #     [-self.data.inst['bl'] / 2,
    #     #      self.data.options.array['ratio'] * self.data.inst['bl'] / 2., 1.]
    #     # ])    Not necessary for get_snr


# -------- Amplitude and abberations methods

# def coupling_efficiency_calculator(self, m: int, alpha: np.ndarray, beta: np.ndarray,
#                                    tip: float = None, tilt: float = None, focus: float = 1.0, iterator: int = 0):
#     if m == 0:
#         if tip or tilt:
#             rho = []
#             a_couple = []
#             for i, wl in enumerate(self.data.inst['wl_bins']):
#                 w_0 = self.data.inst['fundamental_mode_width'][i]
#                 print('fundamental mode width', w_0)
#                 # print('w_0', w_0)
#                 alpha_wl = alpha   #[i, :, :]
#                 beta_wl = beta     #[i, :, :]
#                 # create telescope and fiber curves
#                 z_gauss = np.exp(-((alpha_wl**2+beta_wl**2)/w_0))  # fundamental mode
#                 if i ==0:
#                     fg.single_map(z_gauss,
#                                   extent=[-np.amax(beta_wl), np.amax(beta_wl) , -np.amax(beta_wl), np.amax(beta_wl)],
#                                   labels=['alpha [rad]', 'beta [rad]'], cmap='pink', save=False, save_name='f100_gaussian '+str(iterator))
#                 z_airy = self.airy_function(alpha=alpha_wl, beta=beta_wl, wl=wl,
#                                             tip=tip, tilt=tilt, focus=focus,
#                                             m=i, iterator=iterator)
#                 # print('gauss', z_gauss, z_gauss.shape)
#                 # print('z_airy shape', z_airy.shape)
#                 # calculate meshgrid of alpha, beta
#                 # calculate norms
#                 gauss_norm = np.sqrt(np.trapz(np.trapz(z_gauss**2, beta_wl[:, 0]), alpha_wl[0, :]))
#                 # print('gauss_norms', gauss_norm)
#                 airy_norm = np.sqrt(np.trapz(np.trapz(z_airy**2, beta_wl[:, 0]), alpha_wl[0, :]))
#                 # print('airy_norms', airy_norm)
#                 if i == 0:
#                     fg.single_map(z_airy*z_gauss,
#                                   extent=[-np.amax(beta_wl), np.amax(beta_wl), -np.amax(beta_wl), np.amax(beta_wl)],
#                                   labels=['alpha [rad]', 'beta [rad]'],
#                                   cmap='summer', save=False, save_name='f100_convolution '+str(iterator))
#                 numerator = np.trapz(np.trapz(z_airy*z_gauss, beta_wl[:, 0]), alpha_wl[0, :])
#                 # print('normerator', numerator)
#                 a_couple_wl = numerator/(airy_norm*gauss_norm)
#                 a_couple.append(a_couple_wl)
#                 rho_wl = a_couple_wl**2*100
#                 rho.append(rho_wl)
#             return rho, a_couple
#     else:
#         rho = [1.0]*len(self.data.inst['wl_bins'])
#         return rho
#
# def airy_function(self, alpha: np.ndarray, beta: np.ndarray, m: int, wl: float,
#                   tip: float, tilt: float, focus: float, zoom: float = 1.0, iterator: int = 0):
#     # add Zernicke perturbation to alpha and beta
#     alpha_tip = alpha / focus + tip     # alpha, beta should already be converted to radians
#     beta_tilt = beta / focus + tilt
#     # calculate the Z
#     z = np.pi*self.data.inst['diameter']*np.sqrt(alpha_tip**2+beta_tilt**2) / wl
#     # TODO add Z function zero = small value
#     # calculate airy z for no central obscuration
#     z_airy = (2 * sp.special.jv(1, z)/z)   # TODO check for squared?
#     if m == 0:
#         figure = True
#         if figure:
#             extent_value = 2 * np.arctan(self.data.inst['diameter']/(2 * focus)) #  self.data.inst['hfov_mas'][np.where(self.data.inst['wl_bins'] == wl)][0]
#             print(extent_value)
#             extent = [-extent_value, extent_value, -extent_value, extent_value]
#             fg.single_map(z_airy, extent=extent, cmap='Reds', labels=['alpha [rad]', 'beta [rad]'],  save=False, save_name='f100_airy '+str(iterator))
#     return z_airy
#
# def zernicke_polynomial(self, n: int, m: int, weight):
#     ''' zernicke polynomial function, cartesian representation
#         Taken from Brug (1997)
#         :param n: int
#             order of polynomial
#         :param m: int
#             rotational order of polynomial
#         :param weight :
#         '''
#
#     # polynomials = {'1': x+y^2}
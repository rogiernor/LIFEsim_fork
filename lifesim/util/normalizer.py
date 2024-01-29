# Used for static normalization functi
import numpy as np
import itertools


def transmission_normalizer(transmission_map, normalization_map):
    """
    Normalizer function.
    :param normalization_map: ndarray
        Input normalization map. This fixes all the weird etceteras that have been popping up and hopefully the
        speed
    :param transmission_map: ndarray
        Input transmission map. NB all normalization values are floats, and therefore independent of
        image_size.
    ----------------------------
    :return: normalized_map
        The normalized version of the input map. It switches out the transmission and wavelength to the
        i = wavelength and j = transmission output format.
    """
    normalized_map = np.array([(transmission_map[:, x, :, :]) /
                               (np.amax(normalization_map[x, :, :])) for x in
                               range(len(normalization_map))])  # TODO check if this is an issue with multiple maps
    # normalized_map = np.array([(transmission_map[:, x, :, :] - np.amin(normalization_map[x, :, :])) /
    #                           (np.amax(normalization_map[x, :, :]) - np.amin(normalization_map[x, :, :])) for x in
    #                            range(len(normalization_map))])   # TODO check if this is an issue with multiple maps
    return normalized_map


def onedimensional_map_normalizer(transmission_map, normalization_map=None):
    """
    Normalizer function 2. It sucks! TODO make better!
    :param normalization_map: ndarray
        Input normalization map. This fixes all the weird etceteras that have been popping up and hopefully the
        speed
    :param transmission_map: ndarray
        Input transmission map. NB all normalization values are floats, and therefore independent of
        image_size.
    ----------------------------
    :return: normalized_map
        The normalized version of the input map.
    """

    normalized_map = np.array([(transmission_map[x, :, :] ) /
         (np.amax(normalization_map[x, :, :])) for x in
         range(len(normalization_map))])

    return normalized_map


def baseline_distances(x: np.ndarray):
    n_apertures = len(x)
    norms = []
    for combination in itertools.combinations(range(n_apertures), 2):
        diff = x[combination[0]] - x[combination[1]]
        norms.append(np.sqrt(diff[0] ** 2 + diff[1] ** 2))
    baseline_lengths = norms
    shortest_dist = min(norms)
    return shortest_dist, baseline_lengths

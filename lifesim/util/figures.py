import sys
from typing import Union
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as ptch
import matplotlib.cm as mcm
from collections.abc import Iterable
import matplotlib
matplotlib.rcParams.update({'font.size': 18})

def transmission_figures(maps: list,
                         slab: int,
                         row_col: tuple,
                         extent: list,
                         save=False,
                         save_name='transmission_figure',
                         plot_title='',
                         cmap=mcm.get_cmap('Reds_r')):
    """
    Used to automatically generate figures of the output transmission maps.

    Can also be used to generate the

    Parameters
    ----------
    maps: list
        a list containing the maps you want to visualize
    slab : int
        indicated wavelength slice
    row_col : tuple
        row, column number of plots. Please fit correctly or you get key_error
    extent: list
        defining axis labels NB DOES NOT CHANGE ASPECT
    optional
        save: bool
            defines if figure is saved
        save_name: string
            Defines figure save name in thesis figures file
        plot_title: string
            Defines figure name -> standard is 'output modes'

    :returns: None but figure shows and is possibly saved

    """

    fig, axes = plt.subplots(nrows=row_col[0], ncols=row_col[1], figsize=(float(row_col[1] * 4), float(row_col[0] * 3)))
    # print(axes)
    # print(type(axes))
    print(len(maps))
    for i, mappert in enumerate(maps):
        print('i', i)
        im = axes[i // 2, i % 2].imshow(mappert[slab].T,
                                        vmin=np.amin(mappert[slab]),
                                        vmax=np.amax(mappert[slab]),
                                        extent=extent,
                                        cmap=cmap)
        axes[i // 2, i % 2].set_title('mode ' + str(i + 1))
        fig.colorbar(im, ax=axes[i // 2, i % 2])

        # fig.subplots_adjust(right=0.8)
    fig.suptitle(plot_title)
    plt.tight_layout()
    plt.show()

    if save == True:
        fig.savefig('C:/Users/rogie/Desktop/Afstudeerproject/05. Figures (made by me)/02. Thesis/'
                    + save_name)

def ttn_figures(maps: list,
                    slab: int,
                    extent: list,
                    save=False,
                    save_name='transmission_figure',
                    plot_title=''):
    """
    Used to automatically generate figures of the output transmission maps.

    Can also be used to generate the

    Parameters
    ----------
    maps: list
        a list containing the maps you want to visualize, with dim [n_map, n_bin, alpha, beta]
    slab : int
        indicated wavelength slice
    extent: list
        defining axis labels NB DOES NOT CHANGE ASPECT
    optional
        save: bool
            defines if figure is saved
        save_name: string
            Defines figure save name in thesis figures file
        plot_title: string
            Defines figure name -> standard is 'output modes'

    :returns: None but figure shows and is possibly saved

    """
    n_kernels = len(maps)//2
    print('n_k', n_kernels)
    if n_kernels == 1:
        fig, axes = plt.subplots(nrows=2, ncols=2,
                                 figsize=(float(2 * 3.5), float(2 * 3)))
        kernel = maps[1][slab] - maps [2][slab]
        i = 0
        while i < 3:
            cmap = mcm.get_cmap('Reds_r')
            im = axes[i // 2, i % 2].imshow(maps[i][slab].T,
                                            vmin=np.amin(maps[i][slab]),
                                            vmax=np.amax(maps[i][slab]),
                                            extent=extent,
                                            cmap=cmap)
            fig.colorbar(im, ax=axes[i // 2, i % 2], label='Transmission [-]')
            axes[i // 2, i % 2].set_ylabel('beta [mas]')
            axes[i // 2, i % 2].set_xlabel('alpha [mas]')
            i += 1
        cmap = 'seismic'
        im = axes[1, 1].imshow(kernel.T,
                                        vmin=np.amin(kernel),
                                        vmax=np.amax(kernel),
                                        extent=extent,
                                        cmap=cmap)
        fig.colorbar(im, ax=axes[1, 1], label='Modulation [-]')
        axes[1, 1].set_ylabel('alpha [mas]')
        axes[1, 1].set_xlabel('beta [mas]')
    fig.suptitle(plot_title)
    plt.tight_layout()
    plt.show()

    if save == True:
        fig.savefig('C:/Users/rogie/Desktop/Afstudeerproject/05. Figures (made by me)/02. Thesis/'
                    + save_name)

def single_map(map: any,
               extent=None,
               save=False,
               save_name='2D_figure',
               plot_title='A map',
               intensity='[-]',
               labels=None,
               zoom=False,
               vmin = None,
               vmax = None,
               cmap='seismic'):
    """
    Used to automatically generate figures of the output transmission maps.

    Can also be used to generate other 2d images

    Parameters
    ----------
    maps: dict
        The 2D array you want to visualize. To display it correctly using imshow, the input array is transposed
    slab : int
        indicated wavelength slice
    row_col : tuple
        row, column number of plots. Please fit correctly or you get key_error
    extent: list
        defining axis labels NB DOES NOT CHANGE ASPECT
    optional
        save: bool
            defines if figure is saved
        save_name: string
            Defines figure save name in thesis figures file
        plot_title: string
            Defines figure name -> standard is 'output modes'

    :returns: None but figure shows and is possibly saved

    """

    if cmap == 'Reds':
        cmap = mcm.get_cmap('Reds_r')
    if extent is None:
        extent = [-int(len(map) / 2), int(len(map) / 2), -int(len(map) / 2), int(len(map) / 2)]
    fig, ax = plt.subplots()
    if vmin == None:
        vmin = np.amin(map, (0, 1))
    if vmax == None:
        vmax = np.amax(map, (0, 1))

    im = ax.imshow(np.transpose(map), cmap=cmap, vmin=vmin, vmax=vmax,
                   extent=extent)  # shading = auto (depreciation warning
    if labels is None:
        ax.set_ylabel('y axis [pix]', fontsize=18)
        ax.set_xlabel('x axis [pix]', fontsize=18)
    else:
        ax.set_ylabel(labels[1], fontsize=18)
        ax.set_xlabel(labels[0], fontsize=18)
    if zoom:
        ax.set_xlim([-20, 20])
        ax.set_ylim([-20, 20])

    fig.colorbar(im, ax=ax, label=intensity)
    fig.tight_layout()
    plt.title(plot_title)
    plt.show()

    if save:
        fig.savefig('C:/Users/rogie/Desktop/Afstudeerproject/05. Figures (made by me)/02. Thesis/'
                    + save_name)


def plot_noise_snr(spectrum, labels: list,
                   noises: Union[np.ndarray, list] = None,
                   snr: Union[np.ndarray, list] = None,
                   signals: Union[np.ndarray, list] = None,
                   noise_bottom: float = 10**-2, save = False, save_name: str = 'noise_plot'):
    fig, ax1 = plt.subplots()
    i = 0
    if isinstance(noises, Iterable):
        for n, noise in enumerate(noises):
             # print('found noise' , noise)
            ax1.semilogy(spectrum, noise, label=labels[i])
            i += 1
    if isinstance(snr, Iterable):
        ax2 = ax1.twinx()
        ax2.plot(spectrum, snr, linestyle='--', color='k', label=labels[i])
        ax2.set_ylabel('SNR [-]')
        ax2.legend(loc='upper right')
        i += 1
    if isinstance(signals, Iterable):
        for m, signal in enumerate(signals):
            stylez = ['dashed', 'dashdot']
            print('i', i, labels[i])
            ax1.semilogy(spectrum, signal, linestyle=stylez[m], color='b', label=labels[i])
            i += 1
    ax1.set_xlabel('wavelength bin [m]')
    ax1.set_ylabel('Flux [ph/s/$\mu$m]')
    ax1.legend(loc='upper left')
    ax1.set_ylim(bottom=noise_bottom, top=10**6)
    ax1.grid()
    plt.show()
    if save:
        fig.savefig('C:/Users/rogie/Desktop/Afstudeerproject/05. Figures (made by me)/02. Thesis/'
                    + save_name)


def configuration_render(configuration: list, diameters, dimensions = None, add_circle: bool = False,
                         save = False, save_name = 'config'):
    # print(configuration)
    x = configuration[:, 0]
    # print('x', x)
    y = configuration[:, 1]
    maximum = np.amax(configuration)
    fig_con, ax_con = plt.subplots(figsize=(6, 6))
    for i, posish in enumerate(x):
        p = ptch.Circle((x[i], y[i]), radius=diameters[i]/2, color='silver')
        ax_con.add_patch(p)
        ax_con.annotate('no. '+str(i+1), (x[i]+abs(0.1*x[i]), y[i]-(0.15*y[i])))
    ax_con.set_aspect('equal', 'box')
    ax_con.set_xlim([-1.2*maximum, 1.2*maximum])  # TODO make relative to configuration
    ax_con.set_ylim([-1.2*maximum,  1.2*maximum])
    if dimensions:
        ax_con.set_xlim([-dimensions,  dimensions])  # TODO make relative to configuration
        ax_con.set_ylim([-dimensions,  dimensions])
    if add_circle:
        radius_of_config = max(np.linalg.norm(configuration, axis=-1))
        circle = ptch.Circle((0,0), radius=radius_of_config, ls='--', color='r', fill=False)
        ax_con.add_patch(circle)
        # ax.legend(fontsize=18)
    ax_con.set_xlabel('x [m]', fontsize=18)
    ax_con.set_ylabel('y [m]', fontsize=18)
    # plt.title('configuration of nulling interferometer apertures')
    fig_con.tight_layout()
    plt.grid()
    plt.show()
    if save:
        fig_con.savefig('C:/Users/rogie/Desktop/Afstudeerproject/05. Figures (made by me)/02. Thesis/'
                    + save_name)


def plotitem_arrow(axs, item, plotted, nx, idx, k, osfrac=0.1, verbose=False,
                   baseoffset=0, linestyle="-", label="X", linewidth=0.025,
                   labels=True, projection="polar", rmax=1., addring=False, zorder=1):
    """
    A function that serves as a macro to plot the complex amplitude vectord for CMP

    axs      : The axis objects for the plot
    item     : The complex amplitude to plot
    plotted  : The complex amplitude already plotted (things will get staggered if
                    two identical complex amplitudes are plotted)
    nx       : The number of columns of plots (indexing is 1D if there is only one column)
    i        : The first index of the plot
    j        : The second index of the plot
    k        : The index of the phasor to plot
    osfrac   : The fraction of the amplitude to use as offset
    verbose  : Gives more informations on what happens
    baseoffset : The offset in the start of the vector to use (only special cases)
    linestyle : Used for plotting dashed lines
    label    : The label to use for the legend
    linewidth : Currently not corrected (the scale is very different from the other plotitem function)
    labels   : Whether to include a little label for each vector
    projection : Whether to use a polar or cartesian projection (cartesian not tested)
    rmax      : The maximum norm to for the plot
    """

    offset = osfrac * np.abs(item) * np.exp(1j * (np.angle(item) + np.pi / 2))
    if k == "black":
        thecolor = "k"
    else:
        thecolor = "C" + str(k)

    if verbose: print("initial", item)
    while item + baseoffset in plotted:
        # item += offset
        baseoffset += offset
        if verbose: print("shifting")
    if projection == "polar":
        a0 = np.angle(baseoffset)
        if verbose: print("final, base angle", a0)
        b0 = np.abs(baseoffset)
        if verbose: print("final, base norm", b0)
        a1 = np.angle(item + baseoffset) - np.angle(baseoffset)
        if verbose: print("final, item angle", a1)
        b1 = np.abs(item + baseoffset) - np.abs(baseoffset)
        if verbose: print("final, item norm", b1)
        a2 = np.angle(offset)
        if verbose: print("final, offset angle", a2)
        b2 = np.abs(offset)
        if verbose: print("final, offset norm", b2)
    else:
        a0 = np.real(baseoffset)
        b0 = np.imag(baseoffset)
        a1 = np.real(item)
        b1 = np.imag(item)
        a2 = np.real(offset)
        b2 = np.imag(offset)
    if nx == 1:
        # axs[i,j].scatter(matrix[i*nx+j,k].real, matrix[i*nx+j,k].imag)
        axs[idx].quiver(a0, b0, a1, b1, scale_units='xy', angles='xy', scale=1,
                        color=thecolor, width=linewidth, headlength=2.5, headaxislength=2.2,
                        linestyle=linestyle, label=label, zorder=zorder)
        if labels:
            axs[idx].text(0.95 * a1, 0.9 * b1, str(k))
        axs[idx].set_aspect("equal")
        axs[idx].set_ylim(0, rmax)
        if addring:
            thetas = np.linspace(0, 2 * np.pi, 100)
            axs[idx].plot(thetas, np.ones_like(thetas) * np.abs(item),
                          color=thecolor, zorder=zorder)

    else:

        # axs[idx].scatter(matrix[i*nx+j,k].real, matrix[i*nx+j,k].imag)
        axs[idx].quiver(a0, b0, a1, b1, scale_units='xy', angles='xy', scale=1,
                        color=thecolor, width=linewidth, headlength=2.5, headaxislength=2.2,
                        linestyle=linestyle, label=label, zorder=zorder)
        if labels:
            axs[idx].text(0.95 * a1, 0.9 * b1, str(k))
        axs[idx].set_aspect("equal")
        axs[idx].set_ylim(0, rmax)
        if addring:
            thetas = np.linspace(0, 2 * np.pi, 100)
            axs[idx].plot(thetas, np.ones_like(thetas) * np.abs(item),
                          color=thecolor, zorder=zorder)
    plotted.append(item + baseoffset)
    return plotted


def plot_outputs_smart(matrix=None, inputfield=None, base_preoffset=None, nx=2, ny=None, legendoffset=(-0.2, 0),
                       # (1.6,0.5),
                       verbose=False, osfrac=0.1, plotsize=3, plotspaces=(0.3, 0.4), onlyonelegend=False,
                       labels=False, legend=True, legendsize=8, legendstring="center left", title=None,
                       projection="polar",
                       out_label=None, rmax=None, show=True, onlyoneticklabel=False, labelsize=15,
                       rlabelpos=20, autorm=False, plotter=plotitem_arrow, mainlinewidth=0.04, outputontop=False,
                       thealpha=0.1, color=("black", "silver"), outlabelloc=None, dpi=100):
    """
    Produces a Complex Matrix Plot (CMP) of a combiner matrix. The matrix represents the phasors in each cell of the matrix. In cases where the matrix is designed to take as an input cophased beams of equal amplitude, the plots can also be seen as a representation of the decomposition of the outputs into the contribution of each input.
    returns a fig, and  axs objects.
    matrix   : The matrix to plot.
    inputfield: An input field fed to the combiner. If provided, a plot of M.dot(np.dag(inputfield))
    nx       : The number of columns of plots
    ny       : The number of rows of plots (if provided)
    legendoffset: An offset for the legend (to put it in a free space)
    verbose  :
    osfrac   : The fraction of the amplitude of phasors to use as offset.
    plotsize : The size (matplotlib inches) of each plot
    plotspaces: The space between plots
    labels   : Whether to use the little numbers on each phasor
    legend   : Whether to use a legend for the colors
    legendsize: The size of the legend
    legendstring: The string used by pyplot to use as reference for the legend location
    title    : A title for the whole figure: either the title string, a None object (use default title) or a False boolean (no title).
    projection: if "polar" will use polar plots
    out_label: colored labels for each of the output plot. Requires an array corresponding to each row of the matrix.
    thealpha : Alpha for output label box
    rmax     : The outer limit of the plot (max amplitude)
    show     : Whether to plt.show the figure at the end.
    onlyoneticklabel: Remove tick labels for all but the bottom left plots
    labelsize: Size fot the tick labels (default: 15))
    rlabelpos: The angle at which to put the amplitude tick labels (default: 20)
    autorm   : Automatically remove empty rows True=auto, False=keep all, Boolean array: the rows to remove
    plotter  : A function for plotting fancy arrow vectors
    """
    special = True

    if inputfield is not None:
        initialmatrix = matrix
        matrix = matrix.dot(np.diag(inputfield))
        outvec = initialmatrix.dot(inputfield)
    else:
        initialmatrix = np.zeros((0, 0))
        outvec = None
    if ny is None:
        ny = matrix.shape[0] // nx
    if matrix.shape[0] % nx != 0:
        ny = ny + 1
    ntot = matrix.shape[0]

    if projection == "polar":
        sharex = "none"
        sharey = "all"
        text_coords = "polar"
    else:
        sharex = "all"
        sharey = "all"
        text_coords = "data"
    if rmax is None:
        rmax = np.max(matrix)
    if base_preoffset is None:
        base_preoffset = np.zeros_like(matrix)

    fig, axs = plt.subplots(ny, nx, sharex=sharex, sharey=sharey,
                            gridspec_kw={'hspace': plotspaces[0], 'wspace': plotspaces[1]},
                            figsize=(plotsize * nx, plotsize * matrix.shape[0] // nx + 0.5),
                            subplot_kw=dict(projection=projection), dpi=dpi)

    for idx, theax in enumerate(axs.flatten()):
        if (idx == 0) or (not labels):
            addlabel = False
        else:
            addlabel = True

        # Plotting the output result (black stuff) on the bottom!
        if (outvec is not None) and ((idx) < matrix.shape[0]) and not outputontop:
            plotted = []
            baseoffset = 0
            item = outvec[idx]
            plotted = plotter(axs.flat, item, plotted, nx, idx, "black", verbose=verbose,
                              osfrac=osfrac, baseoffset=baseoffset, linewidth=mainlinewidth,
                              linestyle="-", label="Output " + str(idx), labels=addlabel,
                              projection=projection, rmax=rmax, addring=True, zorder=1)

        plotted = []
        adjust = []
        for k in range(matrix.shape[1]):
            if (idx) < matrix.shape[0]:
                item = matrix[idx, k]  # base_preoffset[idx,k]
                baseoffset = base_preoffset[idx, k]
                if item == 0:
                    continue
                # Here we use plotter, the optional function for plotting vectors
                plotted = plotter(axs.flat, item, plotted, nx, idx, k, verbose=verbose,
                                  osfrac=osfrac, baseoffset=baseoffset, linewidth=mainlinewidth,
                                  linestyle="-", label="Input " + str(k), labels=addlabel,
                                  projection=projection, rmax=rmax)

        plotted2 = []
        adjus2t = []
        # Plotting the dashed lines for the matrix itself
        for k in range(initialmatrix.shape[1]):
            if (idx) < initialmatrix.shape[0]:
                item = initialmatrix[idx, k]
                baseoffset = 0
                if item == 0:
                    continue
                if idx != 0:  # We just don't plot the reference for the bright
                    plotted = plotitem(axs.flat, item, plotted2, nx, idx, k,
                                       osfrac=osfrac, baseoffset=baseoffset,
                                       linestyle="--", label=None, labels=False,
                                       projection=projection, rmax=rmax, linewidth=3, zorder=0)
        # Plotting the output result (black stuff)
        if (outvec is not None) and ((idx) < matrix.shape[0]) and outputontop:
            print("we do plot on top")
            baseoffser = 0
            plotted = []
            item = outvec[idx]
            plotted = plotter(axs.flat, item, plotted, nx, idx, "black", verbose=verbose,
                              osfrac=osfrac, baseoffset=baseoffset, linewidth=mainlinewidth,
                              linestyle="-", label="Output", labels=addlabel,
                              projection=projection, rmax=rmax, addring=True, zorder=4)

        if legend:
            if onlyonelegend:
                if idx == 0:
                    axs.flat[idx].legend(loc=legendstring, prop={'size': legendsize}, bbox_to_anchor=legendoffset)
            else:
                axs.flat[idx].legend(loc=legendstring, prop={'size': legendsize}, bbox_to_anchor=legendoffset)

        if out_label is not None:
            # set_trace()
            if outlabelloc is None:
                outlabelloc = (1.18 * np.pi / 2, 1.32 * rmax)
            for idx, theax in enumerate(axs.flatten()):
                # print(theax)
                if color is None:
                    edgecolor = "C" + str((idx) // 2)
                    facecolor = "C" + str((idx) // 2)
                elif color is False:
                    facecolor = "white"
                    edgecolor = "black"
                else:
                    edgecolor, facecolor = color

                theax.text(outlabelloc[0], outlabelloc[1], str(out_label[idx]), size=15,
                           ha="right", va="top",
                           bbox=dict(boxstyle="round",
                                     # facecolor="none",
                                     alpha=thealpha,
                                     facecolor=facecolor,
                                     edgecolor=edgecolor,
                                     ))

    # eliminating the empty plots
    if autorm is True:
        rowstoremove = np.prod(matrix, axis=1) == 0
    elif autorm is False:
        rowstoremove = np.zeros(matrix.shape[0], dtype=np.bool)
    else:
        rowstoremove = autorm
    # Making a pretty legend for the phase term
    xT = np.arange(0, 2 * np.pi, np.pi / 4)
    # xL=['0',r'$\frac{\pi}{4}$',r'$\frac{\pi}{2}$',r'$\frac{3\pi}{4}$',\
    #        r'$\pi$',r'$\frac{5\pi}{4}$',r'$\frac{3\pi}{2}$',r'$\frac{7\pi}{4}$']
    # xL=['0',r'$\pi/4$',r'$\pi/2$',r'$3\pi/4$',\
    #        r'$\pi$',r'$5\pi/4$',r'$3\pi/2}$',r'$7\pi/4$']
    xT = np.arange(0, 2 * np.pi, np.pi / 2)
    xL = ['0', r'$\pi/2$', r'$\pi$', r'$3\pi/2}$']
    # copying teh ytick labels from the first plot
    # yT = np.linspace(0,0.75*rmax.real,3).round(decimals=1)[1:]
    # yL = [str(yT[b]) for b in np.arange(yT.shape[0])]
    removeticklabels = np.zeros_like(rowstoremove)
    if onlyoneticklabel:
        removeticklabels = np.ones_like(rowstoremove)
        removeticklabels[-nx] = 0
    print("removing the labels:", removeticklabels)
    # print("removing the rows:", rowstoremove)

    for i in np.flip(np.arange(matrix.shape[0])):

        fig.axes[i].set_xticks(xT)
        fig.axes[i].set_xticklabels(xL)
        # fig.axes[i].set_rgrids(yT,yL)
        fig.axes[i].yaxis.set_tick_params(labelbottom=True)

        fig.axes[i].set_rlabel_position(rlabelpos)
        fig.axes[i].tick_params(labelsize=labelsize)
        # print("adding labels",xL)
        if rowstoremove[i]:
            fig.axes[i].remove()
        if removeticklabels[i]:
            fig.axes[i].set_xticklabels([])
            # fig.axes[i].set_yticklabels([])

    if title is not False:
        if title is None:
            title = "The null configurations\n of all the %d outputs" % (matrix.shape[0])
        fig.suptitle(title)
    fig.tight_layout()
    if show:
        plt.show()
    return fig, axs


def print_latex(symbolic, imaginary_unit="j"):
    """
    A convenience function to print a sympy symbolic expression into LaTeX source code.
    """
    prefix = "\\begin{equation}\n  "
    suffix = "\n\end{equation}"
    print(prefix + sp.latex(symbolic, imaginary_unit=imaginary_unit) + suffix, file=sys.stdout)

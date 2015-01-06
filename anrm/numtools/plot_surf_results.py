import matplotlib.pyplot as plt
import matplotlib.gridspec as mgridspec
import matplotlib.ticker as mticker
import matplotlib
import mpl_toolkits.mplot3d as mplot3d
import numpy
import math
import itertools
import multiprocessing
import pysb.integrate
import pickle

"""
    Display the posterior of an MCMC walk on a 3-D surface.

    Parameters
    ----------
    mcmc : bayessb.MCMC
        MCMC object to display.
    dim0, dim1 : indices of parameters to display
    mask : bool/int, optional
        If True (default) the annealing phase of the walk will be discarded
        before plotting. If False, nothing will be discarded and all points will
        be plotted. If an integer, specifies the number of steps to be discarded
        from the beginning of the walk.
    walk : bool, optional
        If True (default) render the walk positions. If False, do not render it.
    rejects : bool, optional
        If True (default) render each rejected position with an 'X'. If False,
        do not render them.
    step : int, optional
        Render every `step`th positions along the walk. Defaults to 1 (render
        all positions). Useful to improve performance with very long walks.
    square_aspect : bool, optional
        If True (default) the X and Y scales of the plot will be equal, allowing
        for direct comparison of moves in the corresponding parameter axes. If
        False the scales will auto-adjust to fit the data tightly, allowing for
        visualization of the full variance along both axes.
    margin : float, optional
        Fraction of the X and Y ranges to add as padding to the surface, beyond
        the range of the points in the walk. Defaults to 0.1. Negative values
        are allowed.
    bounds0, bounds1 : array-like, optional
        Explicit ranges (min, max) for X and Y axes. Specifying either disables
        `square_aspect`.
    zmin, zmax : float, optional
        Max/min height (posterior value) for the sampled surface, and the limits
        for the Z axis of the plot. Any surface points outside this range will
        not be rendered. Defaults to the actual range of posterior values from
        the walk and the sampled surface.
    position_base : array-like, optional
        Vector in log10-parameter space providing values for dimensions *other*
        than dim0/dim1 when calculating the posterior surface (values at
        position dim0 and dim1 will be ignored). Defaults to the median of all
        positions in the walk.
    parallelize : bool, optional
        If True (default), use the multiprocessing module to calculate the
        posterior surface in parallel using all available CPU cores. If False,
        do not parallelize.
    gridsize : int, optional
        Number of points along each axis at which to sample the posterior
        surface. The total number of samples will be `gridsize`**2. Defaults to
        20. Increasing this value will produce a smoother posterior surface at
        the expense of more computational time.
    """

import surface_plot_posteriors as spp

dim0 = spp.dim0
dim1 = spp.dim1
experiment = spp.experiment
bounds0 = spp.b0
bounds1 = spp.b1
gridsize = spp.gs

pos_min, pos_max = [0,0],[0,0]
pos_min[0], pos_max[0] = bounds0
pos_min[1], pos_max[1] = bounds1
p0_vals = numpy.linspace(pos_min[0], pos_max[0], gridsize)
p1_vals = numpy.linspace(pos_min[1], pos_max[1], gridsize)
p0_mesh, p1_mesh = numpy.meshgrid(p0_vals, p1_vals)
    
posterior_mesh = pickle.load(open('%s_dims_%s_%s_surface.pkl' % (experiment, dim0, dim1)))
pmesh_min = numpy.nanmin(posterior_mesh)
pmesh_max = numpy.nanmax(posterior_mesh)

# plot 3-D surface
fig = plt.figure()
ax = fig.gca(projection='3d')
polys = ax.plot_surface(p0_mesh, p1_mesh, posterior_mesh,
                            rstride=1, cstride=1, cmap=matplotlib.cm.jet,
                            linewidth=0.02, alpha=0.2,
                            vmin=pmesh_min, vmax=pmesh_max)
    
    
ax.set_xbound(*bounds0)
ax.set_ybound(*bounds1)
ax.set_xlabel('log10(%s) [dim0]' % spp.pname0)
ax.set_ylabel('log10(%s) [dim1]' % spp.pname1)
ax.set_zlabel('-ln(posterior)')
plt.show()
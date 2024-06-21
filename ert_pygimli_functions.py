import matplotlib.pyplot as plt
import numpy as np
import pygimli as pg
import pygimli.meshtools as mt
from pygimli.physics import ert
import matplotlib as mpl
import pandas as pd
import matplotlib as mpl
import numpy as np
from matplotlib.colors import LogNorm
import os


#########################################
# Data generation

def create_world(start, end, layers, worldMarker=True):
    """Create a world with the given parameters."""
    return mt.createWorld(start=start, end=end, layers=layers, worldMarker=worldMarker)

def add_block_to_world(pos, radius, marker, boundaryMarker, area):
    """Add a block to the world."""
    
    return mt.createCircle(pos=pos, radius=radius, marker=marker,
                            boundaryMarker=boundaryMarker, area=area) 

def add_polygon_to_world(pos, isClosed=True, addNodes=3, interpolate='spline', marker=5):
    """Add a polygon to the world."""

    return mt.createPolygon(pos, isClosed, addNodes, interpolate, marker)

def create_geometry(world, blocks, polygon):
    """Create the geometry with the given world and blocks."""
    return world + blocks + polygon

def create_scheme(elecs, schemeName):
    """Create a scheme with the given electrodes and scheme name."""
    scheme = ert.createData(elecs=elecs, schemeName=schemeName)
    return scheme

def create_mesh(scheme, geom, rhomap, quality):
    """Create a mesh with the given scheme, geometry, resistivity map and quality."""
    # for p in scheme.sensors():
    #     geom.createNode(p)
    #     geom.createNode(p - [0, 0.1])

    for p in scheme.sensorPositions():
        geom.createNode(p)
        geom.createNode(p + pg.RVector3(0, -0.1))
    
    mesh = mt.createMesh(geom, quality=quality)

    # pg.show(mesh, data=rhomap, label=pg.unit('res'), showMesh=True) 

    return mesh

def generate_data(mesh, scheme, rhomap):
    """Generate data with the given mesh, scheme and resistivity map."""

    # data = ert.simulate(mesh, scheme=scheme, res=rhomap, noiseLevel=0,
    #                     noiseAbs=0)#, seed=1337)
    data = ert.simulate(mesh, scheme=scheme, res=rhomap, noiseLevel=1,
                        noiseAbs=1e-6, seed=1337)
    
    # pg.info(np.linalg.norm(data['err']), np.linalg.norm(data['rhoa']))
    # pg.info('Simulated data', data)
    # pg.info('The data contains:', data.dataMap().keys())

    # pg.info('Simulated rhoa (min/max)', min(data['rhoa']), max(data['rhoa']))
    # pg.info('Selected data noise %(min/max)', min(data['err'])*100, max(data['err'])*100)

    return data

def filter_data(data, threshold = 0):
    """Filter the data with the given threshold."""
    data.markInvalid(data('rhoa') < threshold)
    return data

def save_data(data, filename):
    """Export the data to a file."""
    data.save(filename)

def show_data(data):
    """Show the data."""
    ert.show(data)

#########################################
# Inversion

def initialize_ERTMANAGER(path):
    """Initialize the ERTManager with the given path."""
    return ert.ERTManager(path)

def invert_data(data, mesh, manager, limit=1.5, significant=1):
    """Invert the data with the given parameters."""

    # REPRIS DE ALEXIS
    # data['err'] = ert.estimateError(data, absoluteError=0.001, relativeError=0.03)
    # inv = manager.invert(data, paraDX=0.3,  maxCellArea=3, lam=20,verbose=True, paraDepth=0, quality=33.6)

    inv = manager.invert(lam=20, verbose=True)

    # np.testing.assert_approx_equal(manager.inv.chi2(), limit, significant=significant)
    return inv

def show_result_and_fit(manager):
    """Show the result and fit of the manager."""
    manager.showResultAndFit()

def generate_meshPD(manager):
    return pg.Mesh(manager.paraDomain)

#########################################
# Inversion custom fit

def create_custom_inversion_domain(x, y, marker=2):
    """Create an inversion domain with the given parameters."""
    return pg.createGrid(x=x, y=y, marker=marker)

def create_custom_grid(inversionDomain, marker, xbound, ybound):
    """Create a custom grid with the given parameters."""
    return pg.meshtools.appendTriangleBoundary(inversionDomain, marker=marker,
                                               xbound=xbound, ybound=ybound)

def invert_custom_data(manager, grid, data, lam=20, verbose=True):
    """Invert the custom data with the given parameters."""
    return manager.invert(data=data, mesh=grid, lam=lam, verbose=verbose)

#########################################
# Visualization

import matplotlib.gridspec as gridspec
import warnings


def visualize_inversion(model_id, data, manager, mesh, rhomap, rho_min, rho_max):
    # Create figure with manually adjusted subplot positions
    fig = plt.figure(constrained_layout=False, figsize=(14, 5.5))

    # Define positions for each subplot
    pos_left = [0.05, 0.3, 0.3, 0.578]  # [left, bottom, width, height] for the left subplot
    pos_top_right = [0.4, 0.55, 0.55, 0.35]  # [left, bottom, width, height] for the top-right subplot
    pos_bottom_right = [0.4, 0.1, 0.55, 0.35]  # [left, bottom, width, height] for the bottom-right subplot

    # Add subplots with the defined positions
    ax_left = fig.add_axes(pos_left)
    ax_top_right = fig.add_axes(pos_top_right, sharex=ax_left)
    ax_bottom_right = fig.add_axes(pos_bottom_right, sharex=ax_top_right)

    # Plot data
    pg.show(data, ax=ax_left, cMap="jet", logScale=True, cMin=rho_min, cMax=rho_max, colorBar=False)
    ax_left.set_title('Apparent resistivity (Wenner-Schlumberger)')

    # Plot initial model
    pg.show(mesh, rhomap, ax=ax_top_right, hold=True, cMap="jet", logScale=True, colorBar=False, cMin=rho_min, cMax=rho_max)
    ax_top_right.set_xlim([-2.5, 180])
    ax_top_right.set_ylim([-40, 0])
    ax_top_right.set_title('True model')

    # Plot inverse model
    manager.showResult(ax=ax_bottom_right, cMap="jet", logScale=True, cMin=rho_min, cMax=rho_max, colorBar=False)
    ax_bottom_right.set_xlim([-2.5, 180])
    ax_bottom_right.set_ylim([-40, 0])
    ax_bottom_right.set_title('Inversion unstructured grid')

    # Create an axes for the colorbar underneath the left plot
    cbar_ax = fig.add_axes([0.05, 0.18, 0.3, 0.02])  # [left, bottom, width, height] for the colorbar

    # Create the colorbar with logarithmic scale
    cmap = mpl.cm.jet
    norm = mpl.colors.LogNorm(vmin=rho_min, vmax=rho_max)
    cb1 = mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=norm, orientation='horizontal')

    # Add a label to the colorbar
    cb1.set_label('Resistivity [Ohm·m]', rotation=0, labelpad=5)

    # Set centered suptitle
    fig.suptitle('Model ID : ' + model_id + ' - WS, 72 electrodes, 2.5m spacing', x=0.5)

    fig.savefig('figure_results/' + model_id + '.png')
    # plt.close(fig)



def final_plot(model_id, manager, mesh, meshPD, rhomap, inv, synth_data, rho_min, rho_max):


    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize = ((10, 12)))


    pg.show(mesh, rhomap, ax=ax1, hold=True, cMap="jet", logScale=True, orientation="vertical", cMin=rho_min, cMax=rho_max)
    # pg.show(meshPD, ax = ax2, cMap = 'jet', logScale=True, orientation='vertical', cMin=rho_min, cMax=rho_max)
    # pg.show(meshPD, inv, ax=ax2, hold=True, cMap="jet", logScale=True, orientation="vertical", cMin=rho_min, cMax=rho_max)
    manager.showResult(ax=ax2, cMin=rho_min, cMax=rho_max, orientation="vertical")

    # ax1.set_xlim([-50, 50])
    # ax1.set_ylim([-50, 0])
    # ax2.set_xlim([-50, 50])
    # ax2.set_ylim([-50, 0])

    fig.show()


def exportFile(data_geo2x, spacing = 10, name='name.ohm'):
    ndata = len(data_geo2x)
    outputfilename = name
    nelectrod = len(np.unique(np.c_[data_geo2x.iloc[:,1], data_geo2x.iloc[:,2]]))
    
    f = open(outputfilename,'w')
    f.write("{:d}# Number of sensors\n".format(nelectrod))
    f.write("#x\tz\n")
    for i in range(nelectrod):
        f.write("{:f}\t{:f}\n".format(i * spacing,0))
    f.write("{:d}# Number of data\n".format(ndata))    
    f.write("#a\tb\tm\tn\trhoa\n")
    for i in range( len(data_geo2x) ):
        a = int(data_geo2x.iloc[i,1])
        b = int(data_geo2x.iloc[i,2])
        m = int(data_geo2x.iloc[i,3])
        n = int(data_geo2x.iloc[i,4])
        R = 10
        f.write("{}\t{}\t{}\t{}\t{}\n".format(a,b,m,n,R))
    f.close()


def implement_geometry(col):
    rhomap = [[1, col[2]],
             [2, col[3]],
             [3, col[4]]]
    
    world = create_world(start=[-20, 0], end=[201, -80], layers=[-col[0], -col[1]], worldMarker=True)

    # elecs = np.linspace(start= -97.5, stop = 97.5, num=78)
    
    
    return rhomap, world
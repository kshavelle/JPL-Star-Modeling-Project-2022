#This file allows you to run one or all of the componets of the distribution. Instructions on how to run specific parts of the file are in the README file in this directory.
from __future__ import print_function
#Package libraries
from sf3dmodels import Model, Plot_model #Model functions
from sf3dmodels.outflow import OutflowModel    #Model functions
from sf3dmodels import Resolution as Res
import sf3dmodels.utils.units as u #Units
import sf3dmodels.rt as rt  #Writing functions for radiative transfer
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sf3dmodels.grid import Overlap

#Extra libraries
from matplotlib import colors
import numpy as np
import os
import time
import argparse

###########################################

#Initiliaze Constants

t0 = time.time()
dx_grid = 1*u.au
au  = 1.49598e13 #1 au is eqaul to 1.5e13 cm

#Grid parameters
sizex = sizey = sizez = 1000 * u.au #Determines grid size.
#This parameter allows the 2d plane plots to terminate at a specific au value, ie: this is what controls where the sphere terminates! See README and master.inf file for info.
#This determines the size of the grid for the 2d plane plots. Determines where the geometry terminates. This determines the X and Y range of the plot. 
#The largest this value can be is 30000 au. Tested this and at 35000 returns the error zero-size array to reduction operation maximum which has no identity.
Nx = Ny = Nz = 150 #Number of divisions for each axis.
GRID = Model.grid([sizex, sizey, sizez], [Nx, Ny, Nz], rt_code = 'radmc3d') #Enabeling radmc3d by setting rt_code equal to radmc3d. Important for writing input files and running radmc3d mctherm later.
NPoints = GRID.NPoints #Number of nodes in the grid.

# Star parameters for all distributions
#Parameters taken from sf3d models example disc+outflow.
#Parameters taken from sf3d models example disc+outflow.
#Location:/Users/islagarraway/kshavelle-JPL2022/star-forming-regions/examples/disc+outflow File: make_disc.py
#LStar, MStar, TStar, RStar.
MStar = 0.57 * u.MSun #Mass of star.
RStar = 43.0 * u.RSun #Radius of star.
LStar = 0.605e4 * u.LSun #Luminosity of star.
MRate = 1.0e-5 * u.MSun_yr #Mass accretion rate of the star.
TStar = 1.34 * u.Tsun #Tempearature of the star.
Rd = 68.0 * u.au #Centrifugal radius of the star.

# Gas to Dus Ratios
gtd = 100 #Gas to dust ratio.
gtdratio = Model.gastodust(gtd, NPoints)


################################

#Define Disk

def plot_disk(dimension='3D'):

    #Inputs:
        # 3D (string) - creates a 3D disk model
        # 2D (string) - creates a 2D disk model
        
    #Outputs:
        #Rdisc, H0, H0sf_cm

    # Disk Parameters
    #Parameters taken from sf3d models example disc+outflow.
    #Location:/Users/islagarraway/kshavelle-JPL2022/star-forming-regions/examples/disc+outflow File: make_disc.py
    H0 = 0.01 * RStar #Scaleheight at RStar (cm).
    H0sf = 0.03 #Disc scale-height factor (H0 = H0sf * RStar) (au).
    Rdisc = 1.0 * Rd #The radius of the disk at the centrifugal radius of the star (cm).
    Rho0 = Res.Rho0(MRate, Rd, MStar) #Normalization density.
    Rtap = 36 * u.au #Radius where the tapering begins (cm).
    rho_s = 1.0e6 * 1e6 / 5 #From cgs to SI. Density at sonic radius.
    Arho = 15.0 #Density scaling factor.
    #rho_min - The minimnum radius value for the disk.
    #rdisc_max - The maximum radius value for the disk, in this case, set to the centrifugal radius of the star.

    #Make density model for disk
    disk_density = Model.density_Hamburgers(RStar, H0sf, Rd, Rho0, Arho, GRID, rho_min = 1.0,
                                       discFlag = True, Rt = Rtap, rdisc_max = Rdisc)
    #Model.density_Hamburgers
    #Location: /Users/islagarraway/kshavelle-JPL2022/star-forming-regions/sf3dmodels File: Model.py

    #Make .dat file
    prop = {'dens_dust': disk_density.total
            }
            
    #Writing the input files for radmc3d
    #Scripts taken from sf3d models example shells_powerlaw_radmc3d.
    ##Parameters taken from sf3d models example disc+outflow.
    #Location:/Users/islagarraway/kshavelle-JPL2022/star-forming-regions/examples/shells_powerlaw_radmc3d File: make_shellsplaw.py
    #For comparisons between the radmc3d input files and the files made by this file, see master.inf.
    radmc = rt.Radmc3d(GRID, nphot=1000000)
    radmc.submodel(prop, output = 'disc.dat')
    wavelength_intervals = [1e-1,5e2,1e4]
    wavelength_divisions = [20,20]
    radmc.write_radmc3d_control(incl_dust=1, setthreads=8, incl_freefree=0)
    radmc.write_amr_grid()
    radmc.write_dust_density(prop['dens_dust']) #Needs mass density (not number density).
    radmc.write_stars(nstars=1, pos=[[0,0,0]], rstars = [u.RSun], mstars = [u.MSun], flux = [[-u.TSun]],
    #flux --> if negative, radmc assumes the input number as the blackbody temperature of the star
                      lam = wavelength_intervals, nxx = wavelength_divisions)
    radmc.write_wavelength_micron(lam = wavelength_intervals, nxx = wavelength_divisions) #lam --> wavelengths in microns, nxx --> number of divisions in between wavelengths
    
    if dimension == '3D':
        #Plotting density in 3d of disk
        #Plotting scripts taken from example transition_disc in sf3d models.
        #Location: /Users/islagarraway/kshavelle-JPL2022/star-forming-regions/examples/transition_disc/HD135344B_symmetric
        weight = 250. #In Kelvin (k).
        fig = plt.figure(figsize=(8,6))
        #ax = plt.axes(projection='3d')
        lims = (-100,100) #Controls the size of the plot.
        canvas3d = Plot_model.Canvas3d(fig=fig, ax_kw={'xlim': lims, 'ylim': lims, 'zlim': lims, 'azim': 50, 'elev': 10})
        #Plot_model.Canvas3d
        # Location: /Users/islagarraway/kshavelle-JPL2022/star-forming-regions/sf3dmodels File: Plot_model
        ax = canvas3d.ax #generated with fig.add_axes from matplotlib, hence all the matplotlib functions are available
        sp = canvas3d.scatter_random(GRID, disk_density.total, weight, GRID_unit=u.au, power=0.4, NRand=9000, prop_color=disk_density.total/1e6, prop_min=1e2, #function arguments
                                     marker = '.', cmap = 'hot', s = 15, edgecolors = 'none', vmin = None, vmax = None, norm = colors.LogNorm())
        cbar = plt.colorbar(sp)
        cbar.ax.set_ylabel(r'$n_{e^-}$ [cm$^{-3}$]')
        ax.view_init(elev=45, azim=30) #Plot props can be redifined using the usual matplotlib methods
        #ax.view_init(elev=90, azim=30)
        #ax.view_init(elev=0, azim=30)
        ax.set_xlabel('X (AU)')
        ax.set_ylabel('Y (AU)')
        ax.set_zlabel('Z (AU)')
        plt.savefig('disk_dens_scatt_45deg.png', dpi = 200, bbox_inches='tight')
        #plt.savefig('disk_dens_scatt_0deg.png', dpi = 200, bbox_inches='tight')
        #plt.savefig('disk_dens_scatt_90deg.png', dpi = 200, bbox_inches='tight')
        plt.show()
    
    elif dimension == '2D':

        #2d Plot for disk (mid and vertical)
        #Scripts taken from disc+outflow example in sf3d models.
        #Location:/Users/islagarraway/kshavelle-JPL2022/star-forming-regions/examples/disc+outflow File: make_disc.py
        dens_plot = disk_density.total / 1e6
        norm = colors.LogNorm()
        vmin, vmax = np.array([1e15, 5e19]) / 1e6
        norm = colors.LogNorm(vmin=vmin, vmax=vmax)
        Plot_model.plane2D(GRID, dens_plot, axisunit = u.au, cmap = 'ocean_r', plane = {'z': 0*u.au},
                           norm = norm, colorlabel = r'$[\rm cm^{-3}]$', output = 'disk_dens_midplane.png', show = False)
        #Plot_model.plane2D
        #Location:/Users/islagarraway/kshavelle-JPL2022/star-forming-regions/sf3dmodels File:Plot_model
        vmin, vmax = np.array([1e14, 1e18]) / 1e6
        norm = colors.LogNorm(vmin=vmin, vmax=vmax)
        Plot_model.plane2D(GRID, dens_plot, axisunit = u.au, cmap = 'ocean_r', plane = {'y': 0*u.au},
                           norm = norm, colorlabel = r'$[\rm cm^{-3}]$', output = 'disk_dens_vertical.png', show = False)
                       
    return Rdisc, H0, H0sf

#Define Outflow

def plot_outflow():

    #Output:
        #z_max, z_min

    # Outflow paramters
    pos_c = np.array([0*u.au, 0*u.au, 0]) #Origin of coordinates on the axis.
    axis = np.array([0,0,1]) #Long axis direction, of arbitrary length.
    z_min = 1*u.au #Axis lower limit to compute the physical properties.
    z_max = 100*u.au #Axis upper limit to compute the physical properties.

    #Outflow physical properties
    #Outflow parameters described in... Location:/Users/islagarraway/kshavelle-JPL2022/star-forming-regions/sf3dmodels/outflow File:core.py
    #Parameters taken from sf3d model example disc+outflow.
    #Location:/Users/islagarraway/kshavelle-JPL2022/star-forming-regions/examples/disc+outflow File: make_outflow.py
    w0 = 2.0*u.au #Width parameters, parameters to compute the outflow half-width at each z value.
    eps = 0.45 #Cannot find in the documentation.
    w = [w0, eps] #Cannot find in the documentation.
    T0 = 10000. #Parameters to compute the outflow temperature, temperature at r_min.
    qT = -0.5 #Powerlaws for temperature.
    temp = [T0, qT] #parameters to compute the Jet temperature.
    dens0 = 1.e14 #parameters to compute the Jet number density.
    qn = -0.5 #Parameters to compute the outflow number density.
    dens = [dens0, qn] #Parameters to compute the outflow number density.
    ionfrac = [1.0,-0.5] #Parameters to compute the ionized fraction of gas.
    abund = [1.0e-4, 0] #Parameters to compute the molecular abundance.
    gtd = 100 #Gas to dust ratio, same as what was established ealier.
    v0 = 20 * 1e3 #km/s, Gas speed at `z0`. Speed normalization factor.

    #Making density distribution for the outflow
    #Scripts taken from sf3d models example disc+outflow.
    #Location:/Users/islagarraway/kshavelle-JPL2022/star-forming-regions/examples/disc+outflow File: make_outflow.py
    Outf = OutflowModel(pos_c, axis, z_min, z_max, dx_grid) #Initializing Class with grid parameters.
    # OutflowModel
    #Location:/Users/islagarraway/kshavelle-JPL2022/star-forming-regions/sf3dmodels/outflow File: core.py
    Outf.reynolds86(w, dens, ionfrac, temp, v0, abund, gtd) #Invoking the outflow model from Reynolds et al. 1986.
    # reynolds86
    #Location:/Users/islagarraway/kshavelle-JPL2022/star-forming-regions/sf3dmodels/outflow File:core.py

    #Make .dat file for the outflow
    prop = {'dens_dust': Outf.density_ion,
            }
    radmc = rt.Radmc3d(Outf.GRID)
    radmc.submodel(prop, output = 'outflow.dat')#Creating a .dat file in a sub-directory (Subgrids).

    #Outflow defining the outflow density distribution
    outflow_density = prop['dens_dust']

    #Plotting outflow in 3d
    #Plotting scripts taken from example transition_disc in sf3d models.
    #Location: /Users/islagarraway/kshavelle-JPL2022/star-forming-regions/examples/transition_disc/HD135344B_symmetric
    weight = 5000. #Kelvin
    fig = plt.figure(figsize=(8,6))
    #ax = plt.axes(projection='3d')
    lims = (-100,100)
    canvas3d = Plot_model.Canvas3d(fig=fig, ax_kw={'xlim': lims, 'ylim': lims, 'zlim': lims, 'azim': 50, 'elev': 10})
    #Canvas3d
    #Location:/Users/islagarraway/kshavelle-JPL2022/star-forming-regions/sf3dmodels File:Plot_model.py
    ax = canvas3d.ax #generated with fig.add_axes from matplotlib, hence all the matplotlib functions are available
    sp = canvas3d.scatter_random(Outf.GRID, outflow_density/1e6, weight, GRID_unit=u.au, power=0.4, NRand=9000, prop_color=outflow_density/1e6, prop_min=1e2, #function arguments
                                 marker = '.', cmap = 'hot', s = 15, edgecolors = 'none', vmin = None, vmax = None, norm = colors.LogNorm())
    cbar = plt.colorbar(sp)
    cbar.ax.set_ylabel(r'$n_{e^-}$ [cm$^{-3}$]')
    ax.view_init(elev=0, azim=30) #Plot props can be redifined using the usual matplotlib methods.
    #elve = 90 - edge on, elev = 0 - face on.
    #ax.view_init(elev=0, azim=30)
    #ax.view_init(elev=0, azim=30)
    ax.set_xlabel('X (AU)')
    ax.set_ylabel('Y (AU)')
    ax.set_zlabel('Z (AU)')
    plt.savefig('outflow_dens_scatt_0deg.png', dpi = 200, bbox_inches='tight')
    #plt.savefig('outflow_dens_scatt_0deg.png', dpi = 200, bbox_inches='tight')
    #plt.savefig('outflow_dens_scatt_90deg.png', dpi = 200, bbox_inches='tight')
    plt.show()

    #2D Plot from Outflow
    
    #Making 2d plot for outflow (vertical and midplane) - not included in mannual or examples!
    #dens_plot = outflow_density
    #norm = colors.LogNorm()
    #vmin, vmax = np.array([9e15, 5e19]) / 1e6
    #norm = colors.LogNorm(vmin=vmin, vmax=vmax)
    #Plot_model.plane2D(Outf.GRID, dens_plot, axisunit = u.au, cmap = 'ocean_r', plane = {'z': 0*u.au},
                       #norm = norm, colorlabel = r'$[\rm cm^{-3}]$', output = 'outflow_dens_midplane.png', show = False)
    #vmin, vmax = np.array([5e14, 1e18]) / 1e6
    #norm = colors.LogNorm(vmin=vmin, vmax=vmax)
    #Plot_model.plane2D(Outf.GRID, dens_plot, axisunit = u.au, cmap = 'ocean_r', plane = {'y': 0*u.au},
                       #norm = norm, colorlabel = r'$[\rm cm^{-3}]$', output = 'outflow_dens_vertical.png', show = False)
    #Not possible for outflow, not in documentation, not in examples.
    #Returns the error: 'Struct' object has no attribute 'XYZcentres'
    
    return z_max, z_min

def plot_sphere(dimension='3D'):

    #Sphere Parameters
    #Parameters taken from the sf3d model example shells_powerlaw_radmc3d.
    #Location:/Users/islagarraway/kshavelle-JPL2022/star-forming-regions/examples/shells_powerlaw_radmc3d File: make_shellsplaw.py
    #r_max and min specified by Raghvendra in email from 03/09
    r_max = 40000*u.au #Maximum distance to the center (cm).
    r_min = 100*u.au #Minimum distance to the center (cm).
    rho0 = 1e8*1e6 #[part/m3] Number density at r_min.
    r_rho = [r_min, 200*u.au, 300*u.au, r_max] #r frontiers for density.
    q_rho = [-2.0, -1.5, -2.5] #Powerlaws for density.
    #q = -2.0 #Density powerlaw.
    #dens_e = 1.4e5 * 1e6 #Electron number density, from cgs to SI.
    #rho_mean = 1.e5*1e6 #Mean number density.

    #Making the density distribution for the sphere
    #Scripts for sphere taken from #Plotting scripts taken from example shells_powerlaw_radmc3d in sf3d models.
    #Location: /Users/islagarraway/kshavelle-JPL2022/star-forming-regions/examples/shells_powerlaw_radmc3d
    #Sphere density can be described by multiple models.
    #Location: /Users/islagarraway/kshavelle-JPL2022/star-forming-regions/sf3dmodels File: Model.py
    #Model.density_Constant(r_max, GRID, envDens = dens_e) #For a constant sphere.
    sphere_density = Model.density_PowerlawShells(r_rho, q_rho, rho0, GRID, rho_min = 1.0e4) #For a powerlaw sphere.
    #Unable to plot values lower than rho_min = 1.0e3, recive the error: 'a' cannot be empty unless no samples are taken.

    #Make .dat file for the sphere density distribtuion (using Radmc3d)
    prop = {'dens_dust': sphere_density.total,
            }
    radmc = rt.Radmc3d(GRID, nphot=1000000)
    radmc.submodel(prop, output = 'sphere.dat')#Creating a .dat file in a sub-directory (Subgrids).

    #Plotting sphere in 3d
    if dimension == '3D':
        #Plotting scripts taken from example transition_disc in sf3d models.
        #Location: /Users/islagarraway/kshavelle-JPL2022/star-forming-regions/examples/transition_disc/HD135344B_symmetric
        weight = 250. #Kelvin
        fig = plt.figure(figsize=(8,6))
        #ax = plt.axes(projection='3d')
        lims = (-300,300) #Controls the size of the plot.
        canvas3d = Plot_model.Canvas3d(fig=fig, ax_kw={'xlim': lims, 'ylim': lims, 'zlim': lims, 'azim': 50, 'elev': 10})
        #Canvas3d
        #Location:/Users/islagarraway/kshavelle-JPL2022/star-forming-regions/sf3dmodels File:Plot_model.py
        ax = canvas3d.ax #generated with fig.add_axes from matplotlib, hence all the matplotlib functions are available
        sp = canvas3d.scatter_random(GRID, sphere_density.total, weight, GRID_unit=u.au, power=0.4, NRand=10000, prop_color=sphere_density.total/1e6, prop_min=1e2, #function arguments
                                     marker = '.', cmap = 'hot', s = 15, edgecolors = 'none', vmin = None, vmax = None, norm = colors.LogNorm())
        cbar = plt.colorbar(sp)
        cbar.ax.set_ylabel(r'$n_{e^-}$ [cm$^{-3}$]')
        ax.view_init(elev=0, azim=30) #Plot props can be redifined using the usual matplotlib methods
        #elve = 90 - edge on, elev = 0 - face on.
        #ax.view_init(elev=0, azim=30)
        #ax.view_init(elev=90, azim=30)
        ax.set_xlabel('X (AU)')
        ax.set_ylabel('Y (AU)')
        ax.set_zlabel('Z (AU)')
        plt.savefig('sphere_dens_scatt_0deg.png', dpi = 200, bbox_inches='tight')
        #plt.savefig('sphere_dens_scatt_0deg.png', dpi = 200, bbox_inches='tight')
        #plt.savefig('sphere_dens_scatt_90deg.png', dpi = 200, bbox_inches='tight')
        plt.show()

    elif dimension == '2D':
        #Plotting density for sphere in 2d
        #Scripts taken from disc+outflow example in sf3d models.
        #Location:/Users/islagarraway/kshavelle-JPL2022/star-forming-regions/examples/disc+outflow File: make_disc.py
        dens_plot = sphere_density.total
        vmin, vmax = np.array([1e15, 5e19]) / 1e6
        norm = colors.LogNorm(vmin=vmin, vmax=vmax)
        Plot_model.plane2D(GRID, dens_plot, axisunit = u.au, cmap = 'ocean_r', plane = {'z': 0*u.au},
                           norm = norm, colorlabel = r'$[\rm cm^{-3}]$', output = 'sphere_dens_midplane.png', show = False)
        vmin, vmax = np.array([1e14, 1e18]) / 1e6
        norm = colors.LogNorm(vmin=vmin, vmax=vmax)
        Plot_model.plane2D(GRID, dens_plot, axisunit = u.au, cmap = 'ocean_r', plane = {'y': 0*u.au},
                           norm = norm, colorlabel = r'$[\rm cm^{-3}]$', output = 'sphere_dens_vertical.png', show = False)
        #Plot_model.plane2D
        # Location:/Users/islagarraway/kshavelle-JPL2022/star-forming-regions/sf3dmodels File:Plot_model
        Plot_model.plane2D(GRID, dens_plot, axisunit = u.au, cmap = 'ocean_r', plane = {'y': 0*u.au},
                           norm = norm, colorlabel = r'$[\rm cm^{-3}]$', output = 'sphere_dens_vertical.png', show = False)
    return r_max, r_min

def combine_models(disc_file, sphere_file, outflow_file):

    #Inputs:
        # 3D (string) - creates a 3D disk model
        # 2D (string) - creates a 2D disk model
        
    #Outputs:
        #r_max, r_min

    #Combining the disk, sphere, and outflow density distribution .dat files
    #Scripts taken from sf3d models example disc+outflow.
    #Scripts taken from disc+outflow example in sf3d models.
    #Location:/Users/islagarraway/kshavelle-JPL2022/star-forming-regions/examples/disc+outflow File: overlap_models.py
    columns = ['id', 'x', 'y', 'z', 'dens_dust'] #Combining columns (must be the same for all files).
    overlap = Overlap(GRID) #Functions to overlap submodels either from files or from prop objects into a single regular grid.
    #Overlap(GRID)
    #Location:/Users/islagarraway/kshavelle-JPL2022/star-forming-regions/sf3dmodels/grid File:core.py
    #data2merge = ['disc.dat','sphere.dat', 'outflow.dat'] #Merge the three distributions.
    data2merge = [disc_file, sphere_file, outflow_file]
    finalprop = overlap.fromfiles(columns, submodels = data2merge) #Making the combined density distribtuion.

    #Making a combined density distribution
    density = finalprop['dens_dust'] #Combined density distribution.
    weight = 1 #Weighting, allows points to be visible or obscured. The lower the weight, the more points can be seen.
    #Plotting the combined figure
    fig = plt.figure(figsize=(8,6))
    #ax = plt.axes(projection='3d')
    lims = (-200,200) #Controls the size of the plot.
    canvas3d = Plot_model.Canvas3d(fig=fig, ax_kw={'xlim': lims, 'ylim': lims, 'zlim': lims, 'azim': 50, 'elev': 10})
    #Canvas3d
    #Location:/Users/islagarraway/kshavelle-JPL2022/star-forming-regions/sf3dmodels File: Plot_model.py
    ax = canvas3d.ax #generated with fig.add_axes from matplotlib, hence all the matplotlib functions are available
    sp = canvas3d.scatter_random(GRID, density, weight, GRID_unit=u.au, power=0.4, NRand=10000, prop_color=density, prop_min=1e2, #function arguments
                                 marker = '.', cmap = 'hot', s = 15, edgecolors = 'none', vmin = None, vmax = None, norm = colors.LogNorm())
    cbar = plt.colorbar(sp)
    cbar.ax.set_ylabel(r'$n_{e^-}$ [cm$^{-3}$]')
    ax.view_init(elev=45, azim=30) #Plot props can be redifined using the usual matplotlib methods
    #ax.view_init(elev=90, azim=30)
    #ax.view_init(elev=0, azim=30)
    #ax.view_init(elev=45, azim=30)
    ax.set_xlabel('X (AU)')
    ax.set_ylabel('Y (AU)')
    ax.set_zlabel('Z (AU)')
    plt.savefig('sphere+outflow+disk_dens_scatt_45deg.png', dpi = 200, bbox_inches='tight')
    #sphere+outflow+disk_dens_scatt_90deg.png
    #sphere+outflow+disk_dens_scatt_0deg.png
    #sphere+outflow+disk_dens_scatt_45deg.png
    plt.show()
    
    return

def output_params(is_disk, is_outflow, is_sphere, Rdisc=None, H0=None, H0sf=None, z_max=None, z_min=None, r_max=None, r_min=None):

    
    #Writing params into .out file
    with open("physical-parms.out", "w") as external_file:
    
        add_text = "1) DISK"
        print(add_text, file=external_file)
        print(add_text)
        add_text = "A) Geometry"
        print(add_text, file=external_file)
        print(add_text)
        
        if is_disk:
            add_text = "1) DISK"
            print(add_text, file=external_file)
            print(add_text)
            add_text = "A) Geometry"
            print(add_text, file=external_file)
            print(add_text)
            add_text = "Disk outer radius:"
            scientific_notation = "{:.2e}".format(Rdisc)
            print(add_text, scientific_notation, file=external_file)
            print(add_text, scientific_notation)
            add_text = "Scaleheight at RStar:"
            scientific_notation = "{:.2e}".format(H0)
            print(add_text, scientific_notation, file=external_file)
            print(add_text, scientific_notation)
            H0sf_cm = H0sf*au #Disc scale-height factor.
            add_text = "Disk scale-height factor:"
            scientific_notation = "{:.2e}".format(H0sf_cm)
            print(add_text, scientific_notation, file=external_file)
            print(add_text, scientific_notation)
            #Not included in documentation.
            #add_text = "Equatorial Plane Dust Density at disk inner radius:"
            #inner_disk_dens_cm =
            #print(add_text, inner_disk_dens_cm, file=external_file)
            #print(add_text, inner_disk_dens_cm)
            #add_text = "Equatorial Plane Dust Density at disk outer radius:"
            #outer_disk_dens =
            #print(add_text, outer_disk_dens_cm, file=external_file)
            #print(add_text, outer_disk_dens_cm)
        if is_outflow:
            add_text = "2) OUTFLOW"
            print(add_text, file=external_file)
            print(add_text)
            add_text = "A) Geometry"
            print(add_text, file=external_file)
            print(add_text)
            add_text = "Axis lower limit to compute the physical properties:"
            scientific_notation = "{:.2e}".format(z_min)
            print(add_text, scientific_notation, file=external_file)
            print(add_text, scientific_notation)
            add_text = "Axis upper limit to compute the physical properties:"
            scientific_notation = "{:.2e}".format(z_max)
            print(add_text, scientific_notation, file=external_file)
            print(add_text, scientific_notation)
            add_text = "B) Density"
            print(add_text, file=external_file)
            print(add_text)
            #Not included in documentation.
        if is_sphere:
            add_text = "3) SPHERE"
            print(add_text, file=external_file)
            print(add_text)
            add_text = "A) Geometry"
            print(add_text, file=external_file)
            print(add_text)
            add_text = "Outer radius:"
            scientific_notation = "{:.2e}".format(r_max)
            print(add_text, scientific_notation, file=external_file)
            print(add_text, scientific_notation)
            add_text = "Inner radius:"
            scientific_notation = "{:.2e}".format(r_min)
            print(add_text, scientific_notation, file=external_file)
            print(add_text, scientific_notation)
            add_text = "B) Density"
            print(add_text, file=external_file)
            print(add_text)
            #Not included in documentation.
            #add_text = "Density at sphere inner radius:"
            #inner_sphere_dens =
            #print(add_text, inner_sphere_dens, file=external_file)
            #print(add_text, inner_sphere_dens)
            #add_text = "Density at sphere outer radius:"
            #outer_sphere_dens =
            #print(add_text, outer_sphere_dens, file=external_file)
            #print(add_text, outer_sphere_dens)
            #add_text = ""
        external_file.close()
    return
    
if __name__ == "__main__":
#Enabling componets (or all) of the distribution and plot type to be specifed and selected...
    parser = argparse.ArgumentParser()
    parser.add_argument('-model_type',
                            type=str,
                            help='Input desired model type (disk, outflow, sphere, all')
    parser.add_argument('-dimension',
                            type=str,
                            help='2D or 3D output')
    args = parser.parse_args()
    
    if args.model_type == 'disk':
        Rdisc, H0, H0sf = plot_disk(args.dimension)
        output_params(True, False, False, Rdisc=Rdisc, H0=H0, H0sf=H0sf)
    if args.model_type == 'outflow':
        z_max, z_min = plot_outflow()
        output_params(False, True, False, z_max=z_max, z_min=z_min)
    if args.model_type == 'sphere':
        r_max, r_min = plot_sphere(args.dimension)
        output_params(False, False, True, r_max=r_max, r_min=r_min)
    if args.model_type == 'all':
        Rdisc, H0, H0sf = plot_disk(args.dimension)
        z_max, z_min = plot_outflow()
        r_max, r_min = plot_sphere(args.dimension)
        combine_models('disc.dat','sphere.dat', 'outflow.dat')
        output_params(True, True, True, Rdisc, H0, H0sf, z_max, z_min, r_max, r_min)
    
    
    
    
    

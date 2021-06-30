#### IMPORTS ####
from datetime import datetime
import datetime as dt
import numpy as np
from metpy.units import units
import metpy.calc as mpcalc
import matplotlib.pyplot as plt


def addrefcolorbar(ax,fig,im,clevs):
    '''
    This function adds a new colorbar on its own axes for reflectivity

    Inputs: ax, fig are matplotlib axis/figure objects, im is the contourf object,
    and clevs is the contour levels used in the contourf plot

    Outputs: just call it and it'll put the colorbar in the right place

    This code was adapted from Dr. Kim Wood's Community Tools repo
    '''
    axes_bbox = ax.get_position()
    left = axes_bbox.x0 - 0.0075
    bottom = 0.5
    width = 0.005
    height = 0.37
    cax = fig.add_axes([left, bottom, width, height])
    cbar = plt.colorbar(im, cax=cax, ticks=range(5,80,5),
                        orientation='vertical')
    #cbar.ax.xaxis.set_ticks_position('top')

    cbar.ax.tick_params(labelsize=8)
    cbar.set_label('Composite Reflectivity (dBZ)', size=8)  # MODIFY THIS for other fields!!

def addt2mcolorbar(ax,fig,im,clevs):
    '''
    This function adds a new colorbar on its own axes for CAPE

    Inputs: ax, fig are matplotlib axis/figure objects, im is the contourf object,
    and clevs is the contour levels used in the contourf plot

    Outputs: just call it and it'll put the colorbar in the right place

    This code was adapted from Dr. Kim Wood's Community Tools repo
    '''
    axes_bbox = ax.get_position()
    left = axes_bbox.x0 - 0.0075
    bottom = 0.12
    width = 0.005
    height = 0.37
    cax = fig.add_axes([left, bottom, width, height])
    cbar = plt.colorbar(im, cax=cax, ticks=clevs,
                        orientation='vertical')
    #cbar.ax.xaxis.set_ticks_position('top')

    cbar.ax.tick_params(labelsize=8)
    cbar.set_label('2m Temperature (F)', size=8)  # MODIFY THIS for other fields!!

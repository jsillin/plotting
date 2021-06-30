'''
Add colorbars for various maps
'''
import matplotlib.pyplot as plt

########## SOUNDINGMAPS ##########
def addcapecolorbar(ax,fig,im,clevs):
    '''
    This function adds a new colorbar on its own axes for CAPE

    Inputs: ax, fig are matplotlib axis/figure objects, im is the contourf object,
    and clevs is the contour levels used in the contourf plot

    Outputs: just call it and it'll put the colorbar in the right place

    This code was adapted from Dr. Kim Wood's Community Tools repo
    '''
    axes_bbox = ax.get_position()
    left = axes_bbox.x0
    bottom = 0.17
    width = 0.38
    height = 0.01
    cax = fig.add_axes([left, bottom, width, height])
    cbar = plt.colorbar(im, cax=cax, ticks=clevs, orientation='horizontal')
    #cbar.ax.xaxis.set_ticks_position('top')

    cbar.ax.tick_params(labelsize=8)
    cbar.set_label('Surface-Based CAPE (J/kg)', size=8)  # MODIFY THIS for other fields!!

def addwintercapecolorbar(ax,fig,im,clevs):
    '''
    This function adds a new colorbar on its own axes for CAPE

    Inputs: ax, fig are matplotlib axis/figure objects, im is the contourf object,
    and clevs is the contour levels used in the contourf plot

    Outputs: just call it and it'll put the colorbar in the right place

    This code was adapted from Dr. Kim Wood's Community Tools repo
    '''
    axes_bbox = ax.get_position()
    left = axes_bbox.x0
    bottom = 0.16
    width = 0.8
    height = 0.01
    cax = fig.add_axes([left, bottom, width, height])
    cbar = plt.colorbar(im, cax=cax, ticks=clevs, orientation='horizontal')
    #cbar.ax.xaxis.set_ticks_position('top')

    cbar.ax.tick_params(labelsize=8)
    cbar.set_label('Surface-Based CAPE (J/kg)', size=8)  # MODIFY THIS for other fields!!


def addrefcolorbar(ax,fig,im,clevs):
    '''
    This function adds a new colorbar on its own axes for reflectivity

    Inputs: ax, fig are matplotlib axis/figure objects, im is the contourf object,
    and clevs is the contour levels used in the contourf plot

    Outputs: just call it and it'll put the colorbar in the right place

    This code was adapted from Dr. Kim Wood's Community Tools repo
    '''
    axes_bbox = ax.get_position()
    left = axes_bbox.x0 + 0.39
    bottom = 0.17
    width = 0.38
    height = 0.01
    cax = fig.add_axes([left, bottom, width, height])
    cbar = plt.colorbar(im, cax=cax, ticks=clevs, orientation='horizontal')
    #cbar.ax.xaxis.set_ticks_position('top')

    cbar.ax.tick_params(labelsize=8)
    cbar.set_label('Composite Reflectivity (dBZ)', size=8)  # MODIFY THIS for other fields!!

def addwinterrefcolorbar(ax,fig,im,clevs):
    '''
    This function adds a new colorbar on its own axes for reflectivity

    Inputs: ax, fig are matplotlib axis/figure objects, im is the contourf object,
    and clevs is the contour levels used in the contourf plot

    Outputs: just call it and it'll put the colorbar in the right place

    This code was adapted from Dr. Kim Wood's Community Tools repo
    '''
    axes_bbox = ax.get_position()
    left = axes_bbox.x0 + 0.78
    bottom = 0.2
    width = 0.01
    height = 0.15
    cax = fig.add_axes([left, bottom, width, height])
    cbar = plt.colorbar(im, cax=cax, ticks=clevs, orientation='vertical')
    #cbar.ax.xaxis.set_ticks_position('top')

    cbar.ax.tick_params(labelsize=8)
    cbar.set_label('Rain (dBZ)', size=8)  # MODIFY THIS for other fields!!

def addsnowcolorbar(ax,fig,im,clevs):
    '''
    This function adds a new colorbar on its own axes for reflectivity in winter plots

    Inputs: ax, fig are matplotlib axis/figure objects, im is the contourf object,
    and clevs is the contour levels used in the contourf plot

    Outputs: just call it and it'll put the colorbar in the right place

    This code was adapted from Dr. Kim Wood's Community Tools repo
    '''
    axes_bbox = ax.get_position()
    left = axes_bbox.x0  + 0.78
    bottom = 0.68
    width = 0.01
    height = 0.15
    cax = fig.add_axes([left, bottom, width, height])
    cbar = plt.colorbar(im, cax=cax, ticks=clevs, orientation='vertical')
    #cbar.ax.xaxis.set_ticks_position('top')

    cbar.ax.tick_params(labelsize=8)
    cbar.set_label('Snow (dBZ)', size=8)  # MODIFY THIS for other fields!!


def addsleetcolorbar(ax,fig,im,clevs):
    '''
    This function adds a new colorbar on its own axes for reflectivity in winter plots

    Inputs: ax, fig are matplotlib axis/figure objects, im is the contourf object,
    and clevs is the contour levels used in the contourf plot

    Outputs: just call it and it'll put the colorbar in the right place

    This code was adapted from Dr. Kim Wood's Community Tools repo
    '''
    axes_bbox = ax.get_position()
    left = axes_bbox.x0 + 0.78
    bottom = 0.52
    width = 0.01
    height = 0.15
    cax = fig.add_axes([left, bottom, width, height])
    cbar = plt.colorbar(im, cax=cax, ticks=clevs, orientation='vertical')
    #cbar.ax.xaxis.set_ticks_position('top')

    cbar.ax.tick_params(labelsize=8)
    cbar.set_label('Sleet (dBZ)', size=8)  # MODIFY THIS for other fields!!

def addicecolorbar(ax,fig,im,clevs):
    '''
    This function adds a new colorbar on its own axes for reflectivity in winter plots

    Inputs: ax, fig are matplotlib axis/figure objects, im is the contourf object,
    and clevs is the contour levels used in the contourf plot

    Outputs: just call it and it'll put the colorbar in the right place

    This code was adapted from Dr. Kim Wood's Community Tools repo
    '''
    axes_bbox = ax.get_position()
    left = axes_bbox.x0 + 0.78
    bottom = 0.36
    width = 0.01
    height = 0.15
    cax = fig.add_axes([left, bottom, width, height])
    cbar = plt.colorbar(im, cax=cax, ticks=clevs, orientation='vertical')
    #cbar.ax.xaxis.set_ticks_position('top')

    cbar.ax.tick_params(labelsize=8)
    cbar.set_label('Freezing Rain (dBZ)', size=8)  # MODIFY THIS for other fields!!

########## SEVERE DASHBOARD ##########

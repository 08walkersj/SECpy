#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 17:08:28 2021

@author: simon
"""

#%% Imports
import matplotlib.pyplot as plt
from ideas import polar
from matplotlib.widgets import Button
import matplotlib as mpl
import numpy as np
import functools
import warnings
warnings.filterwarnings("ignore")
from SECS_Vis_Plot_Average_Current import average_current_plot
from SECS_Vis_Plot_Prof import profile_plot
from SECS_Vis_Plot_Props import prop_plot
from Create_Vis_File import Create
data_store= './Plotting_Data/SECS_Visualisation.hdf5'
#%% Create Data File
meridian_path= './SECS_Data/results_singularity_mod_clock2.hdf5'
boundaries_path= './SECS_Data/electrojet_boundaries_Clock.hdf5'
Create(boundaries_path, meridian_path, data_store)
#%% Definitions
#Definitions for Button outputs
class Index(object):
    def back(self, event):
        old_mlt= MLT()
        if old_mlt!='':
            global mlt
            if mlt>0:
                mlt=old_mlt-1
            else:
                mlt=23
            for f in polar_fill: f.remove()
            for f in width_plots: f.remove()
            for f in current_plots: f.remove()
            for f in peak_plots: f.remove()
            for f in profiles: f.remove()
            text.remove()
            Width(season, mlt, *Linear_plots[0], polarax)
            Current(season,mlt, *Linear_plots[1])
            Peak(season,mlt, *Linear_plots[2])
            Profile(clock_angle, season, mlt, *Linear_plots[3])
            Conditons(season, clock_angle, mlt)
            plt.pause(0.1)
            plt.show()
    def forward(self, event):
        old_mlt= MLT()
        if old_mlt!='':
            global mlt
            if mlt<23:
                mlt=old_mlt+1
            else:
                mlt=0
            for f in polar_fill: f.remove()
            for f in width_plots: f.remove()
            for f in current_plots: f.remove()
            for f in peak_plots: f.remove()
            for f in profiles: f.remove()
            text.remove()
            Width(season, mlt, *Linear_plots[0], polarax)
            Current(season,mlt, *Linear_plots[1])
            Peak(season,mlt, *Linear_plots[2])
            Profile(clock_angle, season, mlt, *Linear_plots[3])
            Conditons(season, clock_angle, mlt)
            plt.pause(0.1)
            plt.show()
    def summer(self, event):
        if bsummer.color!='0.95':
            bsummer.hovercolor = bsummer.color
            bsummer.color= '0.95'
            bwinter.hovercolor='0.95'
            bwinter.color='0.85'
            ball.hovercolor='0.95'
            ball.color='0.85'
        global season
        season= 'summer'
        for f in polar_fill: f.remove()
        for f in width_plots: f.remove()
        for f in current_plots: f.remove()
        for f in peak_plots: f.remove()
        for f in profiles: f.remove()
        for f in clock_fill: f.remove()
        text.remove()
        Width(season, mlt, *Linear_plots[0], polarax)
        Current(season,mlt, *Linear_plots[1])
        Peak(season,mlt, *Linear_plots[2])
        Profile(clock_angle, season, mlt, *Linear_plots[3])
        polar_plot(clock_angle, season, polarax, clockax)
        Conditons(season, clock_angle, mlt)
        plt.pause(0.1)
        plt.show()
    def winter(self, event):
        if bwinter.color!='0.95':
            bwinter.hovercolor = bwinter.color
            bwinter.color= '0.95'
            bsummer.hovercolor='0.95'
            bsummer.color='0.85'
            ball.hovercolor='0.95'
            ball.color='0.85'
        global season
        season= 'winter'
        for f in polar_fill: f.remove()
        for f in width_plots: f.remove()
        for f in current_plots: f.remove()
        for f in peak_plots: f.remove()
        for f in profiles: f.remove()
        for f in clock_fill: f.remove()
        text.remove()
        Width(season, mlt, *Linear_plots[0], polarax)
        Current(season,mlt, *Linear_plots[1])
        Peak(season,mlt, *Linear_plots[2])
        Profile(clock_angle, season, mlt, *Linear_plots[3])
        polar_plot(clock_angle, season, polarax, clockax)
        Conditons(season, clock_angle, mlt)
        plt.pause(0.1)
        plt.show()
    def All(self, event):
        if ball.color!='0.95':
            ball.hovercolor = ball.color
            ball.color= '0.95'
            bwinter.hovercolor='0.95'
            bwinter.color='0.85'
            bsummer.hovercolor='0.95'
            bsummer.color='0.85'
        global season
        season= ''
        for f in polar_fill: f.remove()
        for f in width_plots: f.remove()
        for f in current_plots: f.remove()
        for f in peak_plots: f.remove()
        for f in profiles: f.remove()
        for f in clock_fill: f.remove()
        text.remove()
        Width(season, mlt, *Linear_plots[0], polarax)
        Current(season,mlt, *Linear_plots[1])
        Peak(season,mlt, *Linear_plots[2])
        Profile(clock_angle, season, mlt, *Linear_plots[3])
        polar_plot(clock_angle, season, polarax, clockax)
        Conditons(season, clock_angle, mlt)
        plt.pause(0.1)
        plt.show()
xs=[]
ys=[]
# Definitions for clickable plots
def onclick(fig, polar_axis, clock_axis, event):
    old_clock= clock()
    old_mlt= MLT()
    global clock_angle
    global mlt
    #If IMF plot is clicked
    if clock_axis.in_axes(event):
        mlt=old_mlt
        ix, iy = event.xdata, event.ydata
        xs.append(ix)
        ys.append(iy)
        ix= int(round(np.rad2deg(ix)/45)*45)*-1
        for f in width_plots: f.remove()
        for f in current_plots: f.remove()
        for f in peak_plots: f.remove()
        for f in profiles: f.remove()
        for f in clock_fill: f.remove()
        for f in polar_plots: f.remove()
        for f in polar_fill: f.remove()
        text.remove()
        if ix<-180:
            ix+=360
        if iy>.3:
            clock_angle= ix
        else:
            clock_angle=''
        polar_plot(clock_angle, season, polar_axis, clock_axis)
        Width(season, mlt, *Linear_plots[0], polar_axis)
        Current(season,mlt, *Linear_plots[1])
        Peak( season,mlt, *Linear_plots[2])
        Profile(clock_angle, season, mlt, *Linear_plots[3])
        Conditons(season, clock_angle, mlt)
        plt.pause(0.1)
        plt.show()
    #If polar plot is clicked
    elif polar_axis.in_axes(event):
        clock_angle=old_clock
        ix, iy = event.xdata, event.ydata
        mlt, mlat= polar_axis.conv_inv(ix, iy)
        for f in polar_fill: f.remove()
        for f in width_plots: f.remove()
        for f in current_plots: f.remove()
        for f in peak_plots: f.remove()
        for f in profiles: f.remove()
        text.remove()
        if mlat>82:
            mlt=''
        else:
            mlt= round(float(mlt))
        Width(season, mlt, *Linear_plots[0], polar_axis)
        Current(season,mlt, *Linear_plots[1])
        Peak(season,mlt, *Linear_plots[2])
        Profile(clock_angle, season, mlt, *Linear_plots[3])
        Conditons(season, clock_angle, mlt)
        plt.pause(0.1)
        plt.show()
def MLT():
    return mlt
def clock():
    return clock_angle
def polar_plot(clock_angle, season, polar_axis, clock_axis):
    global polar_plots, clock_fill, text
    polar_plots= np.concatenate(average_current_plot(clock_angle, season, data_store, polar_axis))
    if clock_angle!='':
        clock_angle= np.deg2rad(clock_angle*-1)
        clock_fill=[clock_axis.fill_between([clock_angle-np.deg2rad(22.5),
                                       clock_angle+np.deg2rad(22.5)], 
                                        [1], color='orange')]
    else:
        clock_fill=[clock_axis.fill_between([i*np.pi/4,(i+1)*np.pi/4], 
                            [1], color='orange') for i in range(8)]
    plt.pause(0.1)
    plt.show()
def Width(season, mlt, axis1, axis2, polar_axis):
    global width_plots, polar_fill
    width_plots=np.concatenate(prop_plot('Width', season, mlt, data_store, axis1, axis2))
    if mlt!='':
        s= mlt-1
        if s<0:
            s+=24
        polar_fill= [polar_axis.fill_between([s, mlt+1], 
                                             [81], color='orange')]
    else:
        polar_fill= [polar_axis.fill_between([i, i+1],
                                              [81], color='orange') for i in range(24)]

def Current(season, mlt, axis1, axis2):
    global current_plots
    current_plots=np.concatenate(prop_plot('Current', season, mlt, data_store, axis1, axis2))
def Peak(season, mlt, axis1, axis2):
    global peak_plots
    peak_plots=np.concatenate(prop_plot('Peak_Value', season, mlt, data_store, axis1, axis2))
def Profile(clock_angle, season, mlt, axis1, axis2):
    global profiles
    profiles= np.concatenate(profile_plot(clock_angle, season, mlt, data_store, axis1, axis2))
def Conditons(season, clock_angle, mlt):
    global text
    if mlt=='':
        mlt= 'All'
    if season=='':
        season='All'
    if clock_angle=='':
        clock_angle='All'
    text= conditions.text(0.07, 0.2, s= f"Season: {season.capitalize()}\nClock Angle: {clock_angle}\nMLT: {mlt}")
def select_clock(angle, clock_axis):
    global clock_fill
    if angle!='all':
        clock_fill=[clock_axis.fill_between([angle-np.deg2rad(22.5),
                                       angle+np.deg2rad(22.5)], 
                                        [1], color='orange')]
    else:
        clock_fill=[clock_axis.fill_between([i*np.pi/4,(i+1)*np.pi/4], 
                            [1], color='orange') for i in range(8)]
    plt.pause(0.1)
    plt.show()
def select_mlt(mlt, polar_axis):
    global polar_fill
    if mlt!='all':
        s= mlt-1
        if s<0:
            s+=24
        polar_fill= [polar_axis.fill_between([s, mlt+1], 
                                             [81], color='orange')]
    else:
        polar_fill= [polar_axis.fill_between([i, i+1],
                                              [81], color='orange') for i in range(24)]
#%% Creating Figure
fig= plt.figure(num='SECS Visualisation Tool')
fig.suptitle('Divergence-Free SECS Visualisation Tool', size=25)
gs= fig.add_gridspec(7, 6, width_ratios=[0.5, 0.5, 0.29, 0.29, 0.29, 0.13], 
                     height_ratios=[0.2, 0.3, 0.5, 0.5,0.1, 0.2, 0.2], 
                     wspace=0, hspace=0)
#%% Create Linear Subplots
global Linear_plots
Linear_plots= [[fig.add_subplot(gs[:2, 0]),fig.add_subplot(gs[2, 0])],
               [fig.add_subplot(gs[3, 0]),fig.add_subplot(gs[4:, 0])],
               [fig.add_subplot(gs[:2,1]),fig.add_subplot(gs[2,1])],
               [fig.add_subplot(gs[3,1]), fig.add_subplot(gs[4:,1])]]
#Setting up broken subplots
d = .015
for i, (ax, ax2) in enumerate(Linear_plots):
    if i%2==0:
        ax.spines['top'].set_visible(False)
        ax.set_xticks([])
        ax.set_xticklabels([])
        ax2.spines['bottom'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax.tick_params(labeltop=False)  # don't put tick labels at the top
    ax2.xaxis.tick_bottom()
    ax.yaxis.tick_left()
    ax2.yaxis.tick_left()
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
Linear_plots[0][0].set_title('Width')
Linear_plots[0][0].set_ylabel('Width ('+u"\N{DEGREE SIGN}"+'Latitude)', y=0, labelpad=15)
Linear_plots[1][0].set_title('Total Current')
Linear_plots[1][0].set_ylabel('Total Current(' +r"$kA$"+')', x=-0.1, y=0, labelpad=15)
Linear_plots[1][1].set_xlabel('IMF BY (' +r"$nT$"+')', labelpad=15)
Linear_plots[2][0].set_title('Peak Value')
Linear_plots[2][0].set_ylabel('Peak Sheet Current Density(' +r"$Am^{-1}$"+')', y=0, labelpad=15)
Linear_plots[3][0].set_title('Sheet Current Profile')
Linear_plots[3][0].set_ylabel('Sheet Current Density(' +r"$Am^{-1}$"+')', y=0, labelpad=15)
for ax in Linear_plots[3]:
    ax.plot([49, 81], [0]*2, color='black')
#%% Create Polar Subplot
global polarax
polarax= polar(fig.add_subplot(gs[:-2,2:-1], projection= 'polar'))
polarax.ax.set_frame_on(False)
polarax.set_title('Evaluated SECS Polar View')
#%% Create IMF Subplot
global clockax
clockax= fig.add_subplot(gs[0, -1], projection= 'polar')
clockax.set_xticklabels([])
clockax.set_yticklabels([])
clockax.set_frame_on(False)
clockax.set_theta_zero_location('N')
for i in range(0,16, 2):
    clockax.plot([(i-1)*np.pi/8, (i-1)*np.pi/8,(i+1)*np.pi/8, (i+1)*np.pi/8], 
                 [0, 1, 1, 0], color='black')
clockax.set_title('Clock Angle')
#%% Create Colorbar
cax= fig.add_subplot(gs[1:, -1])
cbar=fig.colorbar(mpl.cm.ScalarMappable(norm = mpl.colors.Normalize(vmin=-150, vmax=150), cmap=mpl.cm.seismic), cax= cax)
cbar.set_label('Evaluated Br (nT)', size=20)
#%% Create Buttons
global bback, bforward, bsummer, bwinter, ball
buttons= []
for i in range(3):
    for j in range(2):
        buttons.append(fig.add_subplot(gs[5+j, 2+i]))
# Pick axes and set text
bback= Button(buttons[2], u"\u2190")
bback.label.set_fontsize(50)
bforward= Button(buttons[4], u"\u2192")
bforward.label.set_fontsize(50)
bsummer= Button(buttons[1], 'Summer')
bsummer.label.set_fontsize(10)
bwinter= Button(buttons[3], 'Winter')
bwinter.label.set_fontsize(10)
ball= Button(buttons[5], 'All')
ball.label.set_fontsize(10)
ball.hovercolor = '0.85'
ball.color= '0.95'
#Assign functions to buttons
callback= Index()
bback.on_clicked(callback.back)
bforward.on_clicked(callback.forward)
bsummer.on_clicked(callback.summer)
bwinter.on_clicked(callback.winter)
ball.on_clicked(callback.All)
#%% Set Up Clickable Plots
onclick_wrapper=functools.partial(onclick, fig, polarax, clockax)
cid = fig.canvas.mpl_connect('button_press_event', onclick_wrapper)
#%% Conditions Box
from matplotlib import rc
rc('text', usetex=True)
global conditions
conditions= buttons[0]
conditions.set_yticklabels([])
conditions.set_xticklabels([])
conditions.set_xticks([])
conditions.set_yticks([])
conditions.set_xlim(0,1)
conditions.set_ylim(0,1)
conditions.text(x=0.07, y= 0.78, s=r'\underline{Conditions}:', size=10)
conditions.text(x=0.6, y= 0.8, s=r'\underline{Legend}:', size=10)
conditions.text(x=0.6, y= 0.63, s='--- : BY +ve', color= 'red', size=7)
conditions.text(x=0.6, y= 0.46, s='- - : BY -ve', color= 'red', size=7)
conditions.text(x=0.6, y= 0.29, s='--- : High BY', color= 'blue', size=7)
conditions.text(x=0.6, y= 0.12, s='--- : Low BY', color= 'green', size=7)
#%% Initial Set Up
global mlt, season, clock_angle
mlt=''
season=''
clock_angle=''
polar_plot(clock_angle, season, polarax, clockax)
Width(season, mlt, *Linear_plots[0],polarax)
Current(season, mlt, *Linear_plots[1])
Peak(season, mlt, *Linear_plots[2])
Profile(clock_angle, season, mlt, *Linear_plots[3])
Conditons(season, clock_angle, mlt)
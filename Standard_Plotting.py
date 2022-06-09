#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 18:30:23 2022

@author: simon
"""
import numpy as np
def Cubed_Sphere_df(ax_mag, ax_amp, ax_cur, node_grid, eval_grid,
                    poles, MagLon, MagLat, East, South, Up, mag_scale=(1e-5)/8, cur_scale=(1e1)/2,
                    Br_max=5E2, Cur_max=500, Amp_max=6E-1):
    #Magnetic field
    MagI, MagT= poles.Magnetic_Field(East, South, Up)
    Mag= MagI[0]+ MagT[0], MagI[1]+ MagT[1 ], MagI[2]+ MagT[2]
    xi, eta, MagE, MagN= node_grid.projection.vector_cube_projection(*Mag[1:],eval_grid.lon.flatten(), eval_grid.lat.flatten())
    MagQ1=ax_mag.quiver(xi[::8], eta[::8], MagE[::8], MagN[::8], zorder=80,
                        scale=(1e-5)/8, color='Black', alpha=0.7)
    Magqk1= ax_mag.quiverkey(MagQ1, 0.24, 0.175, 400e-9, '', labelpos='E',
                       coordinates='figure', zorder=500)
    
    BR=ax_mag.pcolormesh(eval_grid.xi_mesh, eval_grid.eta_mesh, (Mag[0].reshape(eval_grid.shape))*1E9, vmin=-Br_max, vmax=Br_max,
                                    cmap='bwr', zorder=50, alpha=.5)
    Measured=[None, None]
    xi, eta, Measured[0], Measured[1]= node_grid.projection.vector_cube_projection(East, -South, MagLon, MagLat)
    MagQ2=ax_mag.quiver(xi, eta, *Measured ,zorder=90, 
                         scale = mag_scale, color='red', alpha=0.7)
    
    Magqk2= ax_mag.quiverkey(MagQ2, 0.24, 0.15, 400e-9, r'$400 \ nT$'+'\n'*2, labelpos='E',
                       coordinates='figure', zorder=500)
    #Amplitude
    Amp=ax_amp.pcolormesh(node_grid.xi_mesh, node_grid.eta_mesh, poles.Amplitude(East, South, Up).reshape(node_grid.shape)/node_grid.A , 
                    vmin= -Amp_max, vmax=Amp_max, cmap='seismic', alpha=0.5)
    #Currents
    CI= poles.Currents(East, South, Up)
    CMag=np.sqrt((CI[0]**2) +CI[1]**2)
    xi, eta, CurE, CurN= node_grid.projection.vector_cube_projection(*CI,eval_grid.lon.flatten(), eval_grid.lat.flatten())
    CurQ=ax_cur.quiver(xi[::3], eta[::3], CurE[::3], CurN[::3], zorder=90,scale=cur_scale, color='Black', alpha=0.7)
    Cur=ax_cur.pcolormesh(eval_grid.xi_mesh, eval_grid.eta_mesh, (CMag.reshape(eval_grid.shape))*1e3, vmin=0, vmax= Cur_max, cmap= 'afmhot_r',  zorder=50, alpha=0.5)
    Curqk= ax_cur.quiverkey(CurQ,.75, 0.17, 500e-3, r'$500 \ mAm^{-1}$', labelpos='E',
               coordinates='figure', zorder=500)
    return Magqk1, Magqk2, Curqk, MagQ1, MagQ2, BR, Amp, CurQ, Cur
def Cubed_Sphere_df_amps(ax_df, ax_cf, ax_scat, ax_dir, node_grid, eval_grid,
                         poles, East, South, Up, Amp_data, date, Apex,
                         cur_scale=(1e1)/2, curkey_val=500, Amp_max=6E-1):
    from scipy.interpolate import griddata
    from pandas import Timestamp    
    from astropy.convolution import convolve, Box1DKernel, Gaussian2DKernel
    #Amplitude
    amps_df= poles.Amplitude(East, South, Up).reshape(node_grid.shape)
    Amp_df=ax_df.pcolormesh(node_grid.xi_mesh, node_grid.eta_mesh, amps_df/node_grid.A, 
                    vmin= -Amp_max, vmax=Amp_max, cmap='seismic', alpha=0.5)
    #Currents
    CI= poles.Currents(East, South, Up)
    CI= CI[0].reshape(eval_grid.shape)[::5, ::5], CI[1].reshape(eval_grid.shape)[::5, ::5]
    xi, eta, CurE, CurN= node_grid.projection.vector_cube_projection(*CI,eval_grid.lon[::5, ::5].flatten(), eval_grid.lat[::5, ::5].flatten())
    CurQ_df=ax_df.quiver(xi, eta, CurE, CurN, zorder=90,scale=cur_scale, color='Black', alpha=0.7)
    Curqk= ax_df.quiverkey(CurQ_df, .4, 0.17, curkey_val*1e-3, f'{curkey_val}'+r'$ \ mAm^{-1}$', labelpos='E',
               coordinates='figure', zorder=500)
    #Amplitude
    Amp_data=Amp_data.assign({"mlat": (("vector"),90 - Amp_data['colat'].values)})
    Amp_data=Amp_data.assign({"mlon": (("vector"),Apex.mlt2mlon(Amp_data['mlt'].values, Timestamp(date)))})
    glat, glon, _= Apex.apex2geo(Amp_data['mlat'].values, Amp_data['mlon'].values, 780)
    Amp_data=Amp_data.assign({"glon": (("vector"),glon)})
    Amp_data=Amp_data.assign({"glat": (("vector"),glat)})
    xi, eta=node_grid.projection.geo2cube(glon, glat)
    Amp_data=Amp_data.assign({"xi": (("vector"),xi)})
    Amp_data=Amp_data.assign({"eta": (("vector"),eta)})
    ind= (np.isfinite(Amp_data.xi.values))&(np.isfinite(Amp_data.eta.values))
    amps=griddata((Amp_data.xi.values[ind], Amp_data.eta.values[ind]), Amp_data.Jr.values[ind]*1e-6, 
                  (node_grid.xi, node_grid.eta))*node_grid.A*1e6
    amps=convolve(amps, Gaussian2DKernel(2))
    Amp_cf=ax_cf.pcolormesh(node_grid.xi_mesh, node_grid.eta_mesh, amps/node_grid.A, 
                    vmin= -Amp_max, vmax=Amp_max, cmap='seismic', alpha=0.5)

    amps_= amps.copy()
    #Currents
    GJe, GJn= poles.eval_G_Matrix_J(singularity_limit=50e3, system='curl-free')
    if len(np.where(~np.isfinite(amps))[0]):
        m= np.max(np.where(~np.isfinite(amps)))
        GJe= GJe.reshape(2601, 50, 50)[:, m:, m:].reshape(2601, (50-m)**2)
        GJn= GJn.reshape(2601, 50, 50)[:, m:, m:].reshape(2601, (50-m)**2)
        print('\n'+'>'*20+f' Unable to regrid onto node grid. Grid has been reduced from index: {np.max(np.where(~np.isfinite(amps)))} '+'<'*20+'\n')
        amps=amps[np.max(np.where(~np.isfinite(amps))):, np.max(np.where(~np.isfinite(amps))):]    
    
    
    
    
    Je, Jn= (GJe@amps.flatten()).reshape(eval_grid.shape)[::5, ::5].flatten(), (GJn@amps.flatten()).reshape(eval_grid.shape)[::5, ::5].flatten()
    xi, eta, CurE, CurN= node_grid.projection.vector_cube_projection(Je, Jn, 
                                                                     eval_grid.lon[::5, ::5].flatten(), 
                                                                     eval_grid.lat[::5, ::5].flatten())
    CurQ_cf=ax_cf.quiver(xi, eta, CurE, CurN, zorder=90,scale=cur_scale, color='Black', alpha=0.7)
    
    #Comparison
    ax_scat.relim()
    scat=ax_scat.scatter(poles.Amplitude(East, South, Up).reshape(node_grid.shape)/node_grid.A,amps_/node_grid.A, color='green')
    amps_[amps_<0]=-1
    amps_[amps_>0]=1
    amps_df[amps_df<0]=-1
    amps_df[amps_df>0]=1
    dif= amps_ -amps_df
    dif[dif<0]=-1
    dif[dif>0]=1
    direction=ax_dir.pcolormesh(node_grid.xi_mesh, node_grid.eta_mesh, dif, cmap='bwr', vmin=-1, vmax=1, alpha=1, zorder=-2)
    return Curqk, CurQ_df, CurQ_cf, Amp_df, Amp_cf, scat, direction
def Cubed_Sphere_full_vec(ax_df, ax_dfcf, ax_cf, node_grid, eval_grid,
                         poles, East, South, Up, Amp_data, date, Apex,
                         cur_scale=(1e1)/2, curkey_val=500, Amp_max=6E-1):
    from scipy.interpolate import griddata
    from pandas import Timestamp
    from astropy.convolution import convolve, Box1DKernel, Gaussian2DKernel

    #Amplitude
    amps_df= poles.Amplitude(East, South, Up).reshape(node_grid.shape)
    Amp_df=ax_df.pcolormesh(node_grid.xi_mesh, node_grid.eta_mesh, amps_df/(node_grid.A), 
                    vmin= -Amp_max, vmax=Amp_max, cmap='seismic', alpha=0.5)
    #Currents
    CI= poles.Currents(East, South, Up)
    CI= CI[0].reshape(eval_grid.shape)[::5, ::5], CI[1].reshape(eval_grid.shape)[::5, ::5]
    xi, eta, CurE, CurN= node_grid.projection.vector_cube_projection(*CI,eval_grid.lon[::5, ::5].flatten(), eval_grid.lat[::5, ::5].flatten())
    CurQ_df=ax_df.quiver(xi, eta, CurE, CurN, zorder=90,scale=cur_scale, color='Black', alpha=0.7)
    Curqk= ax_df.quiverkey(CurQ_df, .5, 0.17, curkey_val*1e-3, f'{curkey_val}'+r'$ \ mAm^{-1}$', labelpos='E',
               coordinates='figure', zorder=500)
    #Amplitude
    Amp_data=Amp_data.assign({"mlat": (("vector"),90 - Amp_data['colat'].values)})
    Amp_data=Amp_data.assign({"mlon": (("vector"),Apex.mlt2mlon(Amp_data['mlt'].values, Timestamp(date)))})
    glat, glon, _= Apex.apex2geo(Amp_data['mlat'].values, Amp_data['mlon'].values, 780)
    Amp_data=Amp_data.assign({"glon": (("vector"),glon)})
    Amp_data=Amp_data.assign({"glat": (("vector"),glat)})
    xi, eta=node_grid.projection.geo2cube(glon, glat)
    Amp_data=Amp_data.assign({"xi": (("vector"),xi)})
    Amp_data=Amp_data.assign({"eta": (("vector"),eta)})
    ind= (np.isfinite(Amp_data.xi.values))&(np.isfinite(Amp_data.eta.values))
    amps=griddata((Amp_data.xi.values[ind], Amp_data.eta.values[ind]), Amp_data.Jr.values[ind]*1e-6, 
                  (node_grid.xi, node_grid.eta))*node_grid.A*1e6
    amps=convolve(amps, Gaussian2DKernel(2))
    Amp_cf=ax_cf.pcolormesh(node_grid.xi_mesh, node_grid.eta_mesh, amps/(node_grid.A), 
                    vmin= -Amp_max, vmax=Amp_max, cmap='seismic', alpha=0.5)
    #Currents
    GJe, GJn= poles.eval_G_Matrix_J(singularity_limit=50e3, system='curl-free')    
    Je, Jn= (GJe@amps.flatten()).reshape(eval_grid.shape)[::5, ::5], (GJn@amps.flatten()).reshape(eval_grid.shape)[::5, ::5]
    xi, eta, CurE, CurN= node_grid.projection.vector_cube_projection(Je, Jn, 
                                                                     eval_grid.lon[::5, ::5], 
                                                                     eval_grid.lat[::5, ::5])
    CurQ_cf=ax_cf.quiver(xi, eta, CurE, CurN, zorder=90,scale=cur_scale, color='Black', alpha=0.7)
    
    #Combined Currents
    xi, eta, CurE, CurN= node_grid.projection.vector_cube_projection(Je+CI[0], Jn+CI[1], 
                                                                     eval_grid.lon[::5, ::5].flatten(), 
                                                                     eval_grid.lat[::5, ::5].flatten())
    CurQ_dfcf=ax_dfcf.quiver(xi, eta, CurE, CurN, zorder=90,scale=cur_scale, color='Black', alpha=0.7)
    return Curqk, CurQ_df, CurQ_dfcf, CurQ_cf, Amp_df, Amp_cf
def North_America_Grid(Apex, MagLon, MagLat, λ1= 1e-23, λ2= 1e-21):
    from secsy import cubedsphere as CS
    from SECpy import SECS, RE
    """ Grid"""
    longlim= -160, -40
    latlim= 35, 85
    f1, f2 = Apex.basevectors_qd(61, -102, 0, coords = 'geo')
    qd_north = f2 / np.linalg.norm(f2)
    East, North= qd_north[0], qd_north[1]
    Gridproj= CS.CSprojection((-102, 61), np.rad2deg(np.arctan2(East, North)))
    node_grid= CS.CSgrid(Gridproj, 8000, 6000, 50, 50)
    xi_e  = np.hstack((node_grid.xi_mesh[0]    , node_grid.xi_mesh [0 , - 1] + node_grid.dxi )) - node_grid.dxi /2 
    eta_e = np.hstack((node_grid.eta_mesh[:, 0], node_grid.eta_mesh[-1,   0] + node_grid.deta)) - node_grid.deta/2 
    eval_grid= CS.CSgrid(Gridproj,
                                   node_grid.L + node_grid.Lres, node_grid.W + node_grid.Wres, 
                                   node_grid.Lres, node_grid.Wres, 
                                   edges = (xi_e, eta_e), R = node_grid.R)
    node_lons, node_lats= node_grid.lon.flatten(), node_grid.lat.flatten()
    
    """ SECS Setup """
    depth=500
    poles= SECS(node_lons, node_lats, eval_grid.lon.flatten(), eval_grid.lat.flatten(), 
                mode='image', image_current_radius=RE-depth*1E3)
    
    # Regularisiation
    Le, Ln=node_grid.get_Le_Ln()
    node_f1, node_f2= Apex.basevectors_qd(node_grid.lat.flatten(), node_grid.lon.flatten(), 110, coords='geo')
    e= node_f1/np.linalg.norm(node_f1, axis=0)
    L= np.diag(e[0]).dot(Le) + np.diag(e[1]).dot(Ln)
    G=poles.G_Matrix(MagLon, MagLat)
    GTG= np.dot(G.T, G)
    matrix= GTG + λ1*np.identity(GTG.shape[0]) +λ2*np.dot(L.T, L)/np.max(np.abs(np.dot(L.T, L)))
    Inverse= poles.Inverse_Matrix(matrix, cond=0)
    poles.fitting_matrix['divergence-free']= np.dot(Inverse,G.T)
    return poles, node_grid, eval_grid
def FennoScandia_Grid(Apex, MagLon, MagLat, λ1= 1e-23, λ2= 1e-21, GICs=False, 
                      reg1=1e-22, reg2=0, **SECS_kwargs):
    from secsy import cubedsphere as CS
    from SECpy import SECS, RE
    if not GICs:
        default_kwargs={'mode':'image', 'image_current_radius':RE-500*1E3}
    else:
        default_kwargs={'mode':'default', 'image_current_radius':RE-1E3}
    SECS_kwargs.update({key:default_kwargs[key] for key in default_kwargs.keys() if not key in SECS_kwargs})
    """ Grid"""
    lon_centre= 17.7
    lat_centre= 68.1
    f1, f2 = Apex.basevectors_qd(lat_centre, lon_centre, 0, coords = 'geo')
    qd_north = f2 / np.linalg.norm(f2)
    East, North= qd_north[0], qd_north[1]
    Gridproj= CS.CSprojection((lon_centre, lat_centre), -np.rad2deg(np.arctan2(East, North)))
    node_grid=CS.CSgrid(Gridproj, 3700, 2200, 50., 50.)
    Evalproj= CS.CSprojection((lon_centre-0.23, lat_centre+0.23), -np.rad2deg(np.arctan2(East, North)))
    eval_grid= CS.CSgrid(Evalproj,  3300, 1800, 50., 50.)
    """ SECS Setup """
    poles= SECS(node_grid.lon.flatten(), node_grid.lat.flatten(), eval_grid.lon.flatten(), eval_grid.lat.flatten(), 
                **SECS_kwargs)
    
    # Regularisiation
    Le, Ln=node_grid.get_Le_Ln()
    node_f1, node_f2= Apex.basevectors_qd(node_grid.lat.flatten(), node_grid.lon.flatten(), 110, coords='geo')
    e= node_f1/np.linalg.norm(node_f1, axis=0)
    L= np.diag(e[0]).dot(Le) + np.diag(e[1]).dot(Ln)
    G=poles.G_Matrix(MagLon, MagLat)
    GTG= np.dot(G.T, G)
    if GICs:
        pole_num= len(node_grid.flatten())
        reg1= reg1*np.identity(GTG.shape)
        reg1[:pole_num, :pole_num]= λ1*np.identity(pole_num)
        reg2= reg2*np.identity(GTG.shape)
        reg2[:pole_num, :pole_num]= λ2*np.dot(L.T, L)/np.max(np.abs(np.dot(L.T, L)))
    else:
        reg1= λ1*np.identity(GTG.shape[0])
        reg2= λ2*np.dot(L.T, L)/np.max(np.abs(np.dot(L.T, L)))
    matrix= GTG +reg1 +reg2
    Inverse= poles.Inverse_Matrix(matrix, cond=0)
    poles.fitting_matrix['divergence-free']= np.dot(Inverse,G.T)
    return poles, node_grid, eval_grid
def three_panel_figure(node_grid, Apex, Br_max=5E3, Cur_max=500, Amp_max=6E-1):
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    """ Figure """
    fig=plt.figure()
    gs= fig.add_gridspec(2, 3, height_ratios=[1, .05])
    ax_mag= fig.add_subplot(gs[0,0])
    ax_mag.set_title('Magnetic Field')
    ax_amp=fig.add_subplot(gs[0,1])
    ax_amp.set_title('Amplitude')
    ax_cur=fig.add_subplot(gs[0, 2])
    ax_cur.set_title('Current')
    for ax in (ax_mag, ax_amp, ax_cur):
        ax.set_xlim(node_grid.xi.min(), node_grid.xi.max())
        ax.set_ylim(node_grid.eta.min(), node_grid.eta.max())
        for cl in node_grid.projection.get_projected_coastlines(resolution = '110m'):
            ax.plot(*cl, color='black',zorder=-1, alpha=.25)
        ax.set_aspect('equal')
        ax.set_axis_off()
    DPI=fig.get_dpi()
    fig.set_size_inches(2000.0/float(DPI),1000.0/float(DPI))
    caxes= [fig.add_subplot(gs[1, i]) for i in range(3)]
    mappable_Br=mpl.cm.ScalarMappable(norm = mpl.colors.Normalize(vmin=-Br_max, vmax=Br_max), cmap=mpl.cm.seismic)
    mappable_Amps=mpl.cm.ScalarMappable(norm = mpl.colors.Normalize(vmin=-Amp_max, vmax=Amp_max), cmap=mpl.cm.seismic)
    mappable_Cur=mpl.cm.ScalarMappable(norm = mpl.colors.Normalize(vmin=0, vmax=Cur_max), cmap=mpl.cm.afmhot_r)
    cbar_Br= fig.colorbar(mappable_Br, cax=caxes[0], orientation='horizontal')
    cbar_Br.set_label(r'Radial Magnetic Field ($nT$)')
    cbar_Amps= fig.colorbar(mappable_Amps, cax=caxes[1], orientation='horizontal')
    cbar_Amps.set_label(r'Current Density ($\mu Am^{-2}$)')
    cbar_Cur= fig.colorbar(mappable_Cur, cax=caxes[2], orientation='horizontal')
    cbar_Cur.set_label(r'Current Magnitude ($Akm^{-1}$)')
    mlat, mlon= Apex.geo2apex(node_grid.lat, node_grid.lon, 0)
    clon = [ax.contour(node_grid.xi, node_grid.eta, mlon, levels = np.sort(np.r_[-180:180:10]), 
                       colors = 'black', alpha=.3) for ax in (ax_mag, ax_amp, ax_cur)]
    clat = [ax.contour(node_grid.xi, node_grid.eta, mlat, levels = np.r_[0:90:5], 
                       colors = 'black', alpha=.3) for ax in (ax_mag, ax_amp, ax_cur)]
    cl= ax_amp.clabel(clon[1], inline=True, fontsize=10, fmt = lambda x: str(int(x)) + '$^\circ$')
    cl2=ax_amp.clabel(clat[1], inline=True, fontsize=10, fmt = lambda x: str(int(x)) + '$^\circ$')
    return fig, (ax_mag, ax_amp, ax_cur)
def fourpanel_df_amp_comparison_figure(node_grid, Apex, Amp_max=6E-1):
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    """ Figure """
    fig=plt.figure()
    gs= fig.add_gridspec(3, 3, height_ratios=[.5, .5, .05], hspace=.5)
    ax_df= fig.add_subplot(gs[:-1, 0])
    ax_df.set_title('Divergence-Free')
    ax_cf= fig.add_subplot(gs[:-1, 1])
    ax_cf.set_title('Curl-Free')
    ax_scat= fig.add_subplot(gs[0, -1])
    ax_scat.set_title('DF-CF Comparison')
    ax_scat.set_xlabel('Divergence-Free')
    ax_scat.set_ylabel('Curl-Free')
    ax_dir= fig.add_subplot(gs[1, -1])
    ax_dir.set_title('Direction Comparison')
    for ax in (ax_df, ax_cf, ax_dir):
        ax.set_xlim(node_grid.xi.min(), node_grid.xi.max())
        ax.set_ylim(node_grid.eta.min(), node_grid.eta.max())
        for cl in node_grid.projection.get_projected_coastlines(resolution = '110m'):
            ax.plot(*cl, color='black',zorder=-1, alpha=.25)
        ax.set_aspect('equal')
        ax.set_axis_off()
    DPI=fig.get_dpi()
    fig.set_size_inches(2000.0/float(DPI),1000.0/float(DPI))
    caxes= [fig.add_subplot(gs[-1, :-1]), fig.add_subplot(gs[-1, -1])]
    mappable_Amps=mpl.cm.ScalarMappable(norm = mpl.colors.Normalize(vmin=-Amp_max, vmax=Amp_max), cmap=mpl.cm.seismic)
    cbar_Amps= fig.colorbar(mappable_Amps, cax=caxes[0], orientation='horizontal')
    cbar_Amps.set_label(r'Current Density ($\mu Am^{-2}$)')
    xlim= caxes[-1].get_xlim()
    ylim= caxes[-1].get_ylim()
    caxes[-1].fill_betweenx(list(ylim), 0, 1/3, color='blue')
    caxes[-1].fill_betweenx(list(ylim), 1/3, 2/3, color='white')
    caxes[-1].fill_betweenx(list(ylim), 2/3, 1, color='red')
    caxes[-1].text(1/6, -.3, s='negative', va='center', ha='center')
    caxes[-1].text(3/6, -.3, s='same', va='center', ha='center')
    caxes[-1].text(5/6, -.3, s='positive', va='center', ha='center')
    
    caxes[-1].set_xlim(xlim)
    caxes[-1].set_ylim(ylim)
    caxes[-1].set_xticks([])
    caxes[-1].set_yticks([])
    mlat, mlon= Apex.geo2apex(node_grid.lat, node_grid.lon, 0)
    clon = [ax.contour(node_grid.xi, node_grid.eta, mlon, levels = np.sort(np.r_[-180:180:10]), 
                       colors = 'black', alpha=.3) for ax in (ax_df, ax_cf)]
    clat = [ax.contour(node_grid.xi, node_grid.eta, mlat, levels = np.r_[0:90:5], 
                       colors = 'black', alpha=.3) for ax in (ax_df, ax_cf)]
    cl= ax_df.clabel(clon[0], inline=True, fontsize=10, fmt = lambda x: str(int(x)) + '$^\circ$')
    cl2=ax_df.clabel(clat[0], inline=True, fontsize=10, fmt = lambda x: str(int(x)) + '$^\circ$')
    return fig, (ax_df, ax_cf, ax_scat, ax_dir)
def threepanel_full_vector_figure(node_grid, Apex, Amp_max=6E-1):
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    fig= plt.figure()
    gs= fig.add_gridspec(2, 3, height_ratios=[1, .05])
    ax_df= fig.add_subplot(gs[0, 0])
    ax_df.set_title('Divergence_Free')
    ax_dfcf= fig.add_subplot(gs[0, 1])
    ax_dfcf.set_title('Divergence-Free + Curl-Free')
    ax_cf= fig.add_subplot(gs[0, 2])
    ax_cf.set_title('Curl-Free')
    cax= fig.add_subplot(gs[1, :])
    for ax in (ax_df, ax_dfcf, ax_cf):
        ax.set_xlim(node_grid.xi.min(), node_grid.xi.max())
        ax.set_ylim(node_grid.eta.min(), node_grid.eta.max())
        for cl in node_grid.projection.get_projected_coastlines(resolution = '110m'):
            ax.plot(*cl, color='black',zorder=-1, alpha=.25)
        ax.set_aspect('equal')
        ax.set_axis_off()
    DPI=fig.get_dpi()
    fig.set_size_inches(2000.0/float(DPI),1000.0/float(DPI))
    mappable_Amps=mpl.cm.ScalarMappable(norm = mpl.colors.Normalize(vmin=-Amp_max, vmax=Amp_max), cmap=mpl.cm.seismic)
    cbar_Amps= fig.colorbar(mappable_Amps, cax=cax, orientation='horizontal')
    cbar_Amps.set_label(r'Current Density ($\mu Am^{-2}$)')
    mlat, mlon= Apex.geo2apex(node_grid.lat, node_grid.lon, 0)
    clon = [ax.contour(node_grid.xi, node_grid.eta, mlon, levels = np.sort(np.r_[-180:180:10]), 
                       colors = 'black', alpha=.3) for ax in (ax_df, ax_dfcf, ax_cf)]
    clat = [ax.contour(node_grid.xi, node_grid.eta, mlat, levels = np.r_[0:90:5], 
                       colors = 'black', alpha=.3) for ax in (ax_df, ax_dfcf, ax_cf)]
    cl= ax_df.clabel(clon[0], inline=True, fontsize=10, fmt = lambda x: str(int(x)) + '$^\circ$')
    cl2=ax_df.clabel(clat[0], inline=True, fontsize=10, fmt = lambda x: str(int(x)) + '$^\circ$')
    return fig, (ax_df, ax_dfcf, ax_cf)
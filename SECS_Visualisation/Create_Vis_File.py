#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:21:25 2021

@author: simon
"""
import pandas as pd
import numpy as np
from progressbar import progressbar
import os
from SECpy import ArgumentError
def boundary_create(inpath, outpath, electro_property, in_format='vaex'):
    if in_format.lower()=='pandas':
        properties_pandas(inpath, outpath, electro_property)
    elif in_format.lower()=='vaex':
        properties_vaex(inpath, outpath, electro_property)
    else:
        raise ArgumentError(f'Invalid in format please choose either pandas or vaex. You chose {in_format}')
def properties_vaex(inpath, outpath, electro_property):
    print('\n>>>>>>',electro_property, '<<<<<<\n')
    import vaex as vx
    in_data= vx.open(inpath) # Loading the vaex file
    mlat= np.linspace(49, 81, 50) # Latitude values where the meridian was evaluated
    UL= mlat[-3] # Upper meridian limit
    LL= mlat[2] # Lower meridian limit
    if os.path.isfile(outpath):
        out_data= vx.open(outpath) # Load then output file if it already exists
        columns= np.array(out_data.get_column_names()) # Load output file columns in order to check what has already been done
        print('Creating Back Up')
        if len(out_path.split('.'))>2:
            back_up_path= ''.join(list(out_path[:1]))+'_back_up.'.join(out_path.split('.')[len(out_path.split('.'))-2:])
        else:
            back_up_path= '_back_up.'.join(out_path.split('.'))
        print(back_up_path)
        out_data.export_hdf5(back_up_path)
    else:
        out_data= vx.from_arrays()
        columns=np.array([], dtype='>U3')
    bias_ind= in_data.Circular_Variance_GSM<0.04
    in_data['BT_GSM']= np.sqrt(in_data.BY_GSM_Mean**2 + in_data.BZ_GSM_Mean**2)
    BT_bins= [np.nanquantile(in_data.BT_GSM.values, i) for i in np.linspace(0, 1, 21)]
    if not np.any(np.char.startswith(columns, 'BT_Bins')):
        out_data['BT_Bins']= np.array([f'{l} - {u}' for l, u in zip(BT_bins[:-1], BT_bins[1:])])
    #Seasonal All mlts
    for months, season in progressbar(zip(np.array([[5,6,7], [11, 12, 1]]), np.array(['summer', 'winter'])), max_value=2, prefix='Seasonal |'):
        prefix= electro_property.lower()+'_'+season
        if np.any(np.char.startswith(columns, prefix)):
            print(f'skipped: {prefix}')
            continue
        season_ind= in_data.Date_UTC.dt.month.isin(months)
        pos_ind=(in_data.BY_GSM_Mean!=9.999e3)&(in_data.BY_GSM>=3)
        mean=[]
        error=[]
        BT= []
        for l, u in zip(BT_bins[:-1], BT_bins[1:]):
            in_data.select(pos_ind&(in_data.BT_GSM>=l)&(in_data.BT_GSM<=u)&bias_ind&season_ind)
            for suffix in ['_Eastward1', '_Westward1']:
                mean.append(float(in_data[electro_property+suffix].mean(selection=True)))
                error.append(float(in_data[electro_property+suffix].std(selection=True)/np.sqrt(in_data[electro_property+suffix].count(selection=True))))
            BT.append(float(in_data.median_approx(in_data.BT_GSM, selection=True)))
        mean=np.array(mean)
        error=np.array(error)
        BT=np.array(BT)
        out_data[prefix+'_Eastward_Pos']= mean[::2]
        out_data[prefix+'_Eastward_Error_Pos']= error[::2]
        out_data[prefix+'_Westward_Pos']= mean[1::2]
        out_data[prefix+'_Westward_Error_Pos']= error[1::2]
        out_data[prefix+'_BT_Pos']= BT
        neg_ind=(in_data.BY_GSM_Mean!=9.999e3)&(in_data.BY_GSM<=-3)
        mean=[]
        error=[]
        BT= []
        for l, u in zip(BT_bins[:-1], BT_bins[1:]):
            in_data.select(neg_ind&(in_data.BT_GSM>=l)&(in_data.BT_GSM<=u)&bias_ind&season_ind)
            for suffix in ['_Eastward1', '_Westward1']:
                mean.append(float(in_data[electro_property+suffix].mean(selection=True)))
                error.append(float(in_data[electro_property+suffix].std(selection=True)/np.sqrt(in_data[electro_property+suffix].count(selection=True))))
            BT.append(float(in_data.median_approx(in_data.BT_GSM, selection=True)))
        mean=np.array(mean)
        error=np.array(error)
        BT=np.array(BT)
        out_data[prefix+'_Eastward_Neg']= mean[::2]
        out_data[prefix+'_Eastward_Error_Neg']= error[::2]
        out_data[prefix+'_Westward_Neg']= mean[1::2]
        out_data[prefix+'_Westward_Error_Neg']= error[1::2]
        out_data[prefix+'_BT_Neg']= BT
        out_data.export_hdf5(outpath)
    for mlt in progressbar(np.arange(0, 24), max_value=24, prefix='MLT Seasonal |'):
        mlt1= mlt-1
        mlt2= mlt+1
        if mlt1<0:
            mlt_ind= (in_data.MLT>=23) | (in_data.MLT<=1)
        else:
            mlt_ind= (in_data.MLT>=mlt1) & (in_data.MLT<=mlt2)
        for months, season in progressbar(zip(np.array([[5,6,7], [11, 12, 1]]), np.array(['summer', 'winter'])), max_value=2, prefix='MLT Seasonal |', suffix=f'| mlt: {mlt}'):
            prefix= electro_property.lower()+'_'+season+f'_{mlt}'
            if np.any(np.char.startswith(columns, prefix)):
                print(f'skipped: {prefix}')
                continue
            season_ind= in_data.Date_UTC.dt.month.isin(months)
            pos_ind=(in_data.BY_GSM_Mean!=9.999e3)&(in_data.BY_GSM>=3)
            mean=[]
            error=[]
            BT= []
            for l, u in zip(BT_bins[:-1], BT_bins[1:]):
                in_data.select(pos_ind&(in_data.BT_GSM>=l)&(in_data.BT_GSM<=u)&bias_ind&mlt_ind&season_ind)
                for suffix in ['_Eastward1', '_Westward1']:
                    # print(in_data[electro_property+suffix].count(selection=True))
                    mean.append(float(in_data[electro_property+suffix].mean(selection=True)))
                    error.append(float(in_data[electro_property+suffix].std(selection=True)/np.sqrt(in_data[electro_property+suffix].count(selection=True))))
                BT.append(float(in_data.median_approx(in_data.BT_GSM, selection=True)))
            mean=np.array(mean)
            error=np.array(error)
            BT=np.array(BT)
            out_data[prefix+'_Eastward_Pos']= mean[::2]
            out_data[prefix+'_Eastward_Error_Pos']= error[::2]
            out_data[prefix+'_Westward_Pos']= mean[1::2]
            out_data[prefix+'_Westward_Error_Pos']= error[1::2]
            out_data[prefix+'_BT_Pos']= BT
            neg_ind=(in_data.BY_GSM_Mean!=9.999e3)&(in_data.BY_GSM<=-3)
            mean=[]
            error=[]
            BT= []
            for l, u in zip(BT_bins[:-1], BT_bins[1:]):
                in_data.select(neg_ind&(in_data.BT_GSM>=l)&(in_data.BT_GSM<=u)&bias_ind&mlt_ind&season_ind)
                for suffix in ['_Eastward1', '_Westward1']:
                    # if mlt>15:
                    #     # print(f'{mlt}: ', in_data.evaluate(in_data[electro_property+suffix], selection=True))
                    #     print(f'{mlt}: ',float(in_data[electro_property+suffix].std(selection=True)))
                    # print(in_data[electro_property+suffix].count(selection=True))
                    # print(np.sqrt(in_data[electro_property+suffix].count(selection=True)))
                    mean.append(float(in_data[electro_property+suffix].mean(selection=True)))
                    error.append(float(in_data[electro_property+suffix].std(selection=True)/np.sqrt(in_data[electro_property+suffix].count(selection=True))))
                BT.append(float(in_data.median_approx(in_data.BT_GSM, selection=True)))
            mean=np.array(mean)
            error=np.array(error)
            BT=np.array(BT)
            out_data[prefix+'_Eastward_Neg']= mean[::2]
            out_data[prefix+'_Eastward_Error_Neg']= error[::2]
            out_data[prefix+'_Westward_Neg']= mean[1::2]
            out_data[prefix+'_Westward_Error_Neg']= error[1::2]
            out_data[prefix+'_BT_Neg']= BT
            out_data.export_hdf5(outpath)
    #All Seasons MLTs
    for mlt in progressbar(np.arange(0, 24), max_value=24, prefix='MLT All Seasons |'):
        mlt1= mlt-1
        mlt2= mlt+1
        if mlt1<0:
            mlt_ind= (in_data.MLT>=23) | (in_data.MLT<=1)
        else:
            mlt_ind= (in_data.MLT>=mlt1) & (in_data.MLT<=mlt2)
        prefix= electro_property.lower()+f'_{mlt}'
        if np.any(np.char.startswith(columns, prefix)):
            print(f'skipped: {prefix}')
            continue
        pos_ind=(in_data.BY_GSM_Mean!=9.999e3)&(in_data.BY_GSM>=3)
        mean=[]
        error=[]
        BT= []
        for l, u in zip(BT_bins[:-1], BT_bins[1:]):
            in_data.select(pos_ind&(in_data.BT_GSM>=l)&(in_data.BT_GSM<=u)&bias_ind&mlt_ind)
            for suffix in ['_Eastward1', '_Westward1']:
                mean.append(float(in_data[electro_property+suffix].mean(selection=True)))
                error.append(float(in_data[electro_property+suffix].std(selection=True)/np.sqrt(in_data[electro_property+suffix].count(selection=True))))
            BT.append(float(in_data.median_approx(in_data.BT_GSM, selection=True)))
        mean=np.array(mean)
        error=np.array(error)
        BT=np.array(BT)
        out_data[prefix+'_Eastward_Pos']= mean[::2]
        out_data[prefix+'_Eastward_Error_Pos']= error[::2]
        out_data[prefix+'_Westward_Pos']= mean[1::2]
        out_data[prefix+'_Westward_Error_Pos']= error[1::2]
        out_data[prefix+'_BT_Pos']= BT
        neg_ind=(in_data.BY_GSM_Mean!=9.999e3)&(in_data.BY_GSM<=-3)
        mean=[]
        error=[]
        BT= []
        for l, u in zip(BT_bins[:-1], BT_bins[1:]):
            in_data.select(neg_ind&(in_data.BT_GSM>=l)&(in_data.BT_GSM<=u)&bias_ind&mlt_ind)
            for suffix in ['_Eastward1', '_Westward1']:
                mean.append(float(in_data[electro_property+suffix].mean(selection=True)))
                error.append(float(in_data[electro_property+suffix].std(selection=True)/np.sqrt(in_data[electro_property+suffix].count(selection=True))))
            BT.append(float(in_data.median_approx(in_data.BT_GSM, selection=True)))
        mean=np.array(mean)
        error=np.array(error)
        BT=np.array(BT)
        out_data[prefix+'_Eastward_Neg']= mean[::2]
        out_data[prefix+'_Eastward_Error_Neg']= error[::2]
        out_data[prefix+'_Westward_Neg']= mean[1::2]
        out_data[prefix+'_Westward_Error_Neg']= error[1::2]
        out_data[prefix+'_BT_Neg']= BT
        out_data.export_hdf5(outpath)
    #ALl Seasons All MLTS
    prefix= electro_property.lower()
    if np.any(np.char.startswith(columns, prefix)):
        print(f'skipped: {prefix}')
    else:
        pos_ind=(in_data.BY_GSM_Mean!=9.999e3)&(in_data.BY_GSM>=3)
        mean=[]
        error=[]
        BT= []
        for l, u in zip(BT_bins[:-1], BT_bins[1:]):
            in_data.select(pos_ind&(in_data.BT_GSM>=l)&(in_data.BT_GSM<=u)&bias_ind)
            for suffix in ['_Eastward1', '_Westward1']:
                mean.append(float(in_data[electro_property+suffix].mean(selection=True)))
                error.append(float(in_data[electro_property+suffix].std(selection=True)/np.sqrt(in_data[electro_property+suffix].count(selection=True))))
            BT.append(float(in_data.median_approx(in_data.BT_GSM, selection=True)))
        mean=np.array(mean)
        error=np.array(error)
        BT=np.array(BT)
        out_data[prefix+'_Eastward_Pos']= mean[::2]
        out_data[prefix+'_Eastward_Error_Pos']= error[::2]
        out_data[prefix+'_Westward_Pos']= mean[1::2]
        out_data[prefix+'_Westward_Error_Pos']= error[1::2]
        out_data[prefix+'_BT_Pos']= BT
        neg_ind=(in_data.BY_GSM_Mean!=9.999e3)&(in_data.BY_GSM<=-3)
        mean=[]
        error=[]
        BT= []
        for l, u in zip(BT_bins[:-1], BT_bins[1:]):
            in_data.select(neg_ind&(in_data.BT_GSM>=l)&(in_data.BT_GSM<=u)&bias_ind)
            for suffix in ['_Eastward1', '_Westward1']:
                mean.append(float(in_data[electro_property+suffix].mean(selection=True)))
                error.append(float(in_data[electro_property+suffix].std(selection=True)/np.sqrt(in_data[electro_property+suffix].count(selection=True))))
            BT.append(float(in_data.median_approx(in_data.BT_GSM, selection=True)))
        mean=np.array(mean)
        error=np.array(error)
        BT=np.array(BT)
        out_data[prefix+'_Eastward_Neg']= mean[::2]
        out_data[prefix+'_Eastward_Error_Neg']= error[::2]
        out_data[prefix+'_Westward_Neg']= mean[1::2]
        out_data[prefix+'_Westward_Error_Neg']= error[1::2]
        out_data[prefix+'_BT_Neg']= BT
        out_data.export_hdf5(outpath)
    #All Seasons Clock Angle All MLTs
    for Clock_Angle in progressbar(np.arange(-180, 190, 45), max_value=9, prefix='Clock Angle |'):
        prefix= electro_property.lower()+f'_Clock{Clock_Angle}'
        if abs(Clock_Angle)== 180:
            clock_ind= (in_data.Clock_GSM_Mean>= np.deg2rad(157.5)) | (in_data.Clock_GSM_Mean<=-np.deg2rad(157.5))
        else:
            clock1= Clock_Angle-22.5
            clock2= Clock_Angle+22.5
            clock_ind= (in_data.Clock_GSM_Mean>=np.deg2rad(clock1)) & (in_data.Clock_GSM_Mean<=np.deg2rad(clock2))
        bad_ind= in_data.BY_GSM_Mean!=9.999e3
        mean=[]
        error=[]
        BT= []
        for l, u in zip(BT_bins[:-1], BT_bins[1:]):
            in_data.select(bad_ind&clock_ind&(in_data.BT_GSM>=l)&(in_data.BT_GSM<=u)&bias_ind)
            for suffix in ['_Eastward1', '_Westward1']:
                mean.append(float(in_data[electro_property+suffix].mean(selection=True)))
                error.append(float(in_data[electro_property+suffix].std(selection=True)/np.sqrt(in_data[electro_property+suffix].count(selection=True))))
            BT.append(float(in_data.median_approx(in_data.BT_GSM, selection=True)))
        mean=np.array(mean)
        error=np.array(error)
        BT=np.array(BT)
        out_data[prefix+'_Eastward']= mean[::2]
        out_data[prefix+'_Eastward_Error']= error[::2]
        out_data[prefix+'_Westward']= mean[1::2]
        out_data[prefix+'_Westward_Error']= error[1::2]
        out_data[prefix+'_BT']= BT
        out_data.export_hdf5(outpath)
    #Seasonal Clock Angle All MLTs
    for Clock_Angle in progressbar(np.arange(-180, 190, 45), max_value=9, prefix='Clock Angle Seasonal |'):
        for months, season in progressbar(zip(np.array([[5,6,7], [11, 12, 1]]), np.array(['summer', 'winter'])), max_value=2, 
                                          prefix='Clock Angle Seasonal |', suffix=f'| Clock: {Clock_Angle}'):
            prefix= electro_property.lower()+'_'+season+f'_Clock{Clock_Angle}'
            if np.any(np.char.startswith(columns, prefix)):
                print(f'skipped: {prefix}')
                continue
            season_ind= in_data.Date_UTC.dt.month.isin(months)
            if abs(Clock_Angle)== 180:
                clock_ind= (in_data.Clock_GSM_Mean>= np.deg2rad(157.5)) | (in_data.Clock_GSM_Mean<=-np.deg2rad(157.5))
            else:
                clock1= Clock_Angle-22.5
                clock2= Clock_Angle+22.5
                clock_ind= (in_data.Clock_GSM_Mean>=np.deg2rad(clock1)) & (in_data.Clock_GSM_Mean<=np.deg2rad(clock2))
            bad_ind= in_data.BY_GSM_Mean!=9.999e3
            mean=[]
            error=[]
            BT= []
            for l, u in zip(BT_bins[:-1], BT_bins[1:]):
                in_data.select(bad_ind&clock_ind&(in_data.BT_GSM>=l)&(in_data.BT_GSM<=u)&bias_ind&season_ind)
                for suffix in ['_Eastward1', '_Westward1']:
                    mean.append(float(in_data[electro_property+suffix].mean(selection=True)))
                    error.append(float(in_data[electro_property+suffix].std(selection=True)/np.sqrt(in_data[electro_property+suffix].count(selection=True))))
                BT.append(float(in_data.median_approx(in_data.BT_GSM, selection=True)))
            mean=np.array(mean)
            error=np.array(error)
            BT=np.array(BT)
            out_data[prefix+'_Eastward']= mean[::2]
            out_data[prefix+'_Eastward_Error']= error[::2]
            out_data[prefix+'_Westward']= mean[1::2]
            out_data[prefix+'_Westward_Error']= error[1::2]
            out_data[prefix+'_BT']= BT
            out_data.export_hdf5(outpath)
    #Seasonal Clock Angle MLTs
    for Clock_Angle in progressbar(np.arange(-180, 190, 45), max_value=9, prefix='Clock Angle Seasonal MLTs |'):
        for mlt in progressbar(np.arange(0, 24), max_value=24, prefix='Clock Angle Seasonal MLTs |', suffix=f'| Clock: {Clock_Angle}'):
            mlt1= mlt-1
            mlt2= mlt+1
            if mlt1<0:
                mlt_ind= (in_data.MLT>=23) | (in_data.MLT<=1)
            else:
                mlt_ind= (in_data.MLT>=mlt1) & (in_data.MLT<=mlt2)
            for months, season in progressbar(zip(np.array([[5,6,7], [11, 12, 1]]), np.array(['summer', 'winter'])), max_value=2, 
                                              prefix='Clock Angle Seasonal MLTs |', suffix=f'| Clock: {Clock_Angle}, mlt: {mlt}'):
                prefix= electro_property.lower()+'_'+season+f'_{mlt}'+f'_Clock{Clock_Angle}'
                if np.any(np.char.startswith(columns, prefix)):
                    print(f'skipped: {prefix}')
                    continue
                season_ind= in_data.Date_UTC.dt.month.isin(months)
                if abs(Clock_Angle)== 180:
                    clock_ind= (in_data.Clock_GSM_Mean>= np.deg2rad(157.5)) | (in_data.Clock_GSM_Mean<=-np.deg2rad(157.5))
                else:
                    clock1= Clock_Angle-22.5
                    clock2= Clock_Angle+22.5
                    clock_ind= (in_data.Clock_GSM_Mean>=np.deg2rad(clock1)) & (in_data.Clock_GSM_Mean<=np.deg2rad(clock2))
                bad_ind= in_data.BY_GSM_Mean!=9.999e3
                mean=[]
                error=[]
                BT= []
                for l, u in zip(BT_bins[:-1], BT_bins[1:]):
                    in_data.select(bad_ind&clock_ind&(in_data.BT_GSM>=l)&(in_data.BT_GSM<=u)&bias_ind&season_ind&mlt_ind)
                    for suffix in ['_Eastward1', '_Westward1']:
                        mean.append(float(in_data[electro_property+suffix].mean(selection=True)))
                        error.append(float(in_data[electro_property+suffix].std(selection=True)/np.sqrt(in_data[electro_property+suffix].count(selection=True))))
                    BT.append(float(in_data.median_approx(in_data.BT_GSM, selection=True)))
                mean=np.array(mean)
                error=np.array(error)
                BT=np.array(BT)
                out_data[prefix+'_Eastward']= mean[::2]
                out_data[prefix+'_Eastward_Error']= error[::2]
                out_data[prefix+'_Westward']= mean[1::2]
                out_data[prefix+'_Westward_Error']= error[1::2]
                out_data[prefix+'_BT']= BT
                out_data.export_hdf5(outpath)
    #All Seasons Clock Angle MLTs
    for Clock_Angle in progressbar(np.arange(-180, 190, 45), max_value=9, prefix='Clock Angle MLTs |'):
        for mlt in progressbar(np.arange(0, 24), max_value=24, prefix='Clock Angle MLTs |', suffix= f'| Clock: {Clock_Angle}'):
            mlt1= mlt-1
            mlt2= mlt+1
            if mlt1<0:
                mlt_ind= (in_data.MLT>=23) | (in_data.MLT<=1)
            else:
                mlt_ind= (in_data.MLT>=mlt1) & (in_data.MLT<=mlt2)
            prefix= electro_property.lower()+f'_{mlt}'+f'_Clock{Clock_Angle}'
            if np.any(np.char.startswith(columns, prefix)):
                print(f'skipped: {prefix}')
                continue
            if abs(Clock_Angle)== 180:
                clock_ind= (in_data.Clock_GSM_Mean>= np.deg2rad(157.5)) | (in_data.Clock_GSM_Mean<=-np.deg2rad(157.5))
            else:
                clock1= Clock_Angle-22.5
                clock2= Clock_Angle+22.5
                clock_ind= (in_data.Clock_GSM_Mean>=np.deg2rad(clock1)) & (in_data.Clock_GSM_Mean<=np.deg2rad(clock2))
            bad_ind= in_data.BY_GSM_Mean!=9.999e3
            mean=[]
            error=[]
            BT= []
            for l, u in zip(BT_bins[:-1], BT_bins[1:]):
                in_data.select(bad_ind&clock_ind&(in_data.BT_GSM>=l)&(in_data.BT_GSM<=u)&bias_ind&mlt_ind)
                for suffix in ['_Eastward1', '_Westward1']:
                    mean.append(float(in_data[electro_property+suffix].mean(selection=True)))
                    error.append(float(in_data[electro_property+suffix].std(selection=True)/np.sqrt(in_data[electro_property+suffix].count(selection=True))))
                BT.append(float(in_data.median_approx(in_data.BT_GSM, selection=True)))
            mean=np.array(mean)
            error=np.array(error)
            BT=np.array(BT)
            out_data[prefix+'_Eastward']= mean[::2]
            out_data[prefix+'_Eastward_Error']= error[::2]
            out_data[prefix+'_Westward']= mean[1::2]
            out_data[prefix+'_Westward_Error']= error[1::2]
            out_data[prefix+'_BT']= BT
            out_data.export_hdf5(outpath)
    # Dipole Tilt All Clocks All MLTs
    for tilt in progressbar(np.arange(-30, 35, 10), max_value=7, prefix='Dipole Tilt All Clock Angles All MLTs'):
        prefix= electro_property.lower()+f'_Tilt{tilt}'
        if np.any(np.char.startswith(columns, prefix)):
            print(f'skipped: {prefix}')
            continue
        tilt1= tilt-5
        tilt2= tilt+5
        tilt_ind= (in_data.Dipole_Tilt>=tilt1) & (in_data.Dipole_Tilt<=tilt2)
        pos_ind=(in_data.BY_GSM_Mean!=9.999e3)&(in_data.BY_GSM>=3)
        mean=[]
        error=[]
        BT= []
        for l, u in zip(BT_bins[:-1], BT_bins[1:]):
            in_data.select(pos_ind&(in_data.BT_GSM>=l)&(in_data.BT_GSM<=u)&bias_ind&tilt_ind)
            for suffix in ['_Eastward1', '_Westward1']:
                mean.append(float(in_data[electro_property+suffix].mean(selection=True)))
                error.append(float(in_data[electro_property+suffix].std(selection=True)/np.sqrt(in_data[electro_property+suffix].count(selection=True))))
            BT.append(float(in_data.median_approx(in_data.BT_GSM, selection=True)))
        mean=np.array(mean)
        error=np.array(error)
        BT=np.array(BT)
        out_data[prefix+'_Eastward_Pos']= mean[::2]
        out_data[prefix+'_Eastward_Error_Pos']= error[::2]
        out_data[prefix+'_Westward_Pos']= mean[1::2]
        out_data[prefix+'_Westward_Error_Pos']= error[1::2]
        out_data[prefix+'_BT_Pos']= BT
        neg_ind=(in_data.BY_GSM_Mean!=9.999e3)&(in_data.BY_GSM<=-3)
        mean=[]
        error=[]
        BT= []
        for l, u in zip(BT_bins[:-1], BT_bins[1:]):
            in_data.select(neg_ind&(in_data.BT_GSM>=l)&(in_data.BT_GSM<=u)&bias_ind&tilt_ind)
            for suffix in ['_Eastward1', '_Westward1']:
                mean.append(float(in_data[electro_property+suffix].mean(selection=True)))
                error.append(float(in_data[electro_property+suffix].std(selection=True)/np.sqrt(in_data[electro_property+suffix].count(selection=True))))
            BT.append(float(in_data.median_approx(in_data.BT_GSM, selection=True)))
        mean=np.array(mean)
        error=np.array(error)
        BT=np.array(BT)
        out_data[prefix+'_Eastward_Neg']= mean[::2]
        out_data[prefix+'_Eastward_Error_Neg']= error[::2]
        out_data[prefix+'_Westward_Neg']= mean[1::2]
        out_data[prefix+'_Westward_Error_Neg']= error[1::2]
        out_data[prefix+'_BT_Neg']= BT
        out_data.export_hdf5(outpath)
    # Dipole Tilt All Clocks MLTs
    for tilt in progressbar(np.arange(-30, 35, 10), max_value=7, prefix='Dipole Tilt All Clock Angles MLTs |'):
        for mlt in progressbar(np.arange(0, 24), max_value=24, prefix='Dipole Tilt All Clock Angles MLTs |', suffix= f'| Tilt: {tilt}'):
            mlt1= mlt-1
            mlt2= mlt+1
            if mlt1<0:
                mlt_ind= (in_data.MLT>=23) | (in_data.MLT<=1)
            else:
                mlt_ind= (in_data.MLT>=mlt1) & (in_data.MLT<=mlt2)
            prefix= electro_property.lower()+f'_{mlt}'+f'_Tilt{tilt}'
            if np.any(np.char.startswith(columns, prefix)):
                print(f'skipped: {prefix}')
                continue
            tilt1= tilt-5
            tilt2= tilt+5
            tilt_ind= (in_data.Dipole_Tilt>=tilt1) & (in_data.Dipole_Tilt<=tilt2)
            pos_ind=(in_data.BY_GSM_Mean!=9.999e3)&(in_data.BY_GSM>=3)
            mean=[]
            error=[]
            BT= []
            for l, u in zip(BT_bins[:-1], BT_bins[1:]):
                in_data.select(pos_ind&(in_data.BT_GSM>=l)&(in_data.BT_GSM<=u)&bias_ind&tilt_ind&mlt_ind)
                for suffix in ['_Eastward1', '_Westward1']:
                    mean.append(float(in_data[electro_property+suffix].mean(selection=True)))
                    error.append(float(in_data[electro_property+suffix].std(selection=True)/np.sqrt(in_data[electro_property+suffix].count(selection=True))))
                BT.append(float(in_data.median_approx(in_data.BT_GSM, selection=True)))
            mean=np.array(mean)
            error=np.array(error)
            BT=np.array(BT)
            out_data[prefix+'_Eastward_Pos']= mean[::2]
            out_data[prefix+'_Eastward_Error_Pos']= error[::2]
            out_data[prefix+'_Westward_Pos']= mean[1::2]
            out_data[prefix+'_Westward_Error_Pos']= error[1::2]
            out_data[prefix+'_BT_Pos']= BT
            neg_ind=(in_data.BY_GSM_Mean!=9.999e3)&(in_data.BY_GSM<=-3)
            mean=[]
            error=[]
            BT= []
            for l, u in zip(BT_bins[:-1], BT_bins[1:]):
                in_data.select(neg_ind&(in_data.BT_GSM>=l)&(in_data.BT_GSM<=u)&bias_ind&tilt_ind&mlt_ind)
                for suffix in ['_Eastward1', '_Westward1']:
                    mean.append(float(in_data[electro_property+suffix].mean(selection=True)))
                    error.append(float(in_data[electro_property+suffix].std(selection=True)/np.sqrt(in_data[electro_property+suffix].count(selection=True))))
                BT.append(float(in_data.median_approx(in_data.BT_GSM, selection=True)))
            mean=np.array(mean)
            error=np.array(error)
            BT=np.array(BT)
            out_data[prefix+'_Eastward_Neg']= mean[::2]
            out_data[prefix+'_Eastward_Error_Neg']= error[::2]
            out_data[prefix+'_Westward_Neg']= mean[1::2]
            out_data[prefix+'_Westward_Error_Neg']= error[1::2]
            out_data[prefix+'_BT_Neg']= BT
            out_data.export_hdf5(outpath)
    # Dipole Tilt Clocks MLTs
    for tilt in progressbar(np.arange(-30, 35, 10), max_value=7, prefix='Dipole Tilt Clock Angles MLTs |'):
        for Clock_Angle in progressbar(np.arange(-180, 190, 45), max_value=9, prefix='Dipole Tilt Clock Angles MLTs |',
                                       suffix=f'| Tilt: {tilt}'):
            for mlt in progressbar(np.arange(0, 24), max_value=24, prefix='Dipole Tilt Clock Angles MLTs |', 
                                   suffix= f'| Tilt: {tilt}, Clock: {Clock_Angle}'):
                prefix= electro_property.lower()+f'_{mlt}'+f'_Clock{Clock_Angle}'+f'_Tilt{tilt}'
                if np.any(np.char.startswith(columns, prefix)):
                    print(f'skipped: {prefix}')
                    continue
                mlt1= mlt-1
                mlt2= mlt+1
                if mlt1<0:
                    mlt_ind= (in_data.MLT>=23) | (in_data.MLT<=1)
                else:
                    mlt_ind= (in_data.MLT>=mlt1) & (in_data.MLT<=mlt2)
                if abs(Clock_Angle)== 180:
                    clock_ind= (in_data.Clock_GSM_Mean>= np.deg2rad(157.5)) | (in_data.Clock_GSM_Mean<=-np.deg2rad(157.5))
                else:
                    clock1= Clock_Angle-22.5
                    clock2= Clock_Angle+22.5
                    clock_ind= (in_data.Clock_GSM_Mean>=np.deg2rad(clock1)) & (in_data.Clock_GSM_Mean<=np.deg2rad(clock2))
                tilt1= tilt-5
                tilt2= tilt+5
                tilt_ind= (in_data.Dipole_Tilt>=tilt1) & (in_data.Dipole_Tilt<=tilt2)
                bad_ind= in_data.BY_GSM_Mean!=9.999e3
                mean=[]
                error=[]
                BT= []
                for l, u in zip(BT_bins[:-1], BT_bins[1:]):
                    in_data.select(bad_ind&(in_data.BT_GSM>=l)&(in_data.BT_GSM<=u)&bias_ind&tilt_ind&mlt_ind&clock_ind)
                    for suffix in ['_Eastward1', '_Westward1']:
                        mean.append(float(in_data[electro_property+suffix].mean(selection=True)))
                        error.append(float(in_data[electro_property+suffix].std(selection=True)/np.sqrt(in_data[electro_property+suffix].count(selection=True))))
                    BT.append(float(in_data.median_approx(in_data.BT_GSM, selection=True)))
                mean=np.array(mean)
                error=np.array(error)
                BT=np.array(BT)
                out_data[prefix+'_Eastward']= mean[::2]
                out_data[prefix+'_Eastward_Error']= error[::2]
                out_data[prefix+'_Westward']= mean[1::2]
                out_data[prefix+'_Westward_Error']= error[1::2]
                out_data[prefix+'_BT']= BT
                out_data.export_hdf5(outpath)
    # Dipole Tilt Clocks All MLTs
    for tilt in progressbar(np.arange(-30, 35, 10), max_value=7, prefix='Dipole Tilt Clock Angles |'):
        for Clock_Angle in progressbar(np.arange(-180, 190, 45), max_value=9, prefix='Dipole Tilt Clock Angles |',
                                       suffix=f'| Tilt: {tilt}'):
            prefix= electro_property.lower()+f'_Clock{Clock_Angle}'+f'_Tilt{tilt}'
            if np.any(np.char.startswith(columns, prefix)):
                print(f'skipped: {prefix}')
                continue
            mlt1= mlt-1
            mlt2= mlt+1
            if mlt1<0:
                mlt_ind= (in_data.MLT>=23) | (in_data.MLT<=1)
            else:
                mlt_ind= (in_data.MLT>=mlt1) & (in_data.MLT<=mlt2)
            if abs(Clock_Angle)== 180:
                clock_ind= (in_data.Clock_GSM_Mean>= np.deg2rad(157.5)) | (in_data.Clock_GSM_Mean<=-np.deg2rad(157.5))
            else:
                clock1= Clock_Angle-22.5
                clock2= Clock_Angle+22.5
                clock_ind= (in_data.Clock_GSM_Mean>=np.deg2rad(clock1)) & (in_data.Clock_GSM_Mean<=np.deg2rad(clock2))
            tilt1= tilt-5
            tilt2= tilt+5
            tilt_ind= (in_data.Dipole_Tilt>=tilt1) & (in_data.Dipole_Tilt<=tilt2)
            bad_ind= in_data.BY_GSM_Mean!=9.999e3
            mean=[]
            error=[]
            BT= []
            for l, u in zip(BT_bins[:-1], BT_bins[1:]):
                in_data.select(bad_ind&(in_data.BT_GSM>=l)&(in_data.BT_GSM<=u)&bias_ind&tilt_ind&clock_ind)
                for suffix in ['_Eastward1', '_Westward1']:
                    mean.append(float(in_data[electro_property+suffix].mean(selection=True)))
                    error.append(float(in_data[electro_property+suffix].std(selection=True)/np.sqrt(in_data[electro_property+suffix].count(selection=True))))
                BT.append(float(in_data.median_approx(in_data.BT_GSM, selection=True)))
            mean=np.array(mean)
            error=np.array(error)
            BT=np.array(BT)
            out_data[prefix+'_Eastward']= mean[::2]
            out_data[prefix+'_Eastward_Error']= error[::2]
            out_data[prefix+'_Westward']= mean[1::2]
            out_data[prefix+'_Westward_Error']= error[1::2]
            out_data[prefix+'_BT']= BT
            out_data.export_hdf5(outpath)
def properties_pandas(inpath, outpath, electro_property):
    print('\n',electro_property, '\n')
    mlat= np.linspace(49, 81, 50)
    UL= mlat[-3]
    LL= mlat[2]
    # Loading the pandas store to check what already exists in the file (avoiding repeats)
    store=pd.HDFStore(inpath, mode='r')
    keys=[]
    if os.path.isfile(outpath):
        out_store= pd.HDFStore(outpath, mode='r')
        keys= out_store.keys()
        out_store.close()
    keys=np.array(keys, dtype='object')
    #Seasonal All mlts
    for months, season in progressbar(zip(np.array([[5,6,7], [11, 12, 1]]), np.array(['summer', 'winter'])), max_value=2):
        key= electro_property+'_'+season
        if np.any(keys=='/'+key):
            continue
        tmp_df= pd.DataFrame()
        df= store.select(key='main', where="BY_GSM!=9.999e3 & Peak_Value_Westward1!=9.999e3 & Poleward_Boundary_Westward1<UL & Equatorward_Boundary_Westward1>LL",
                         columns=['Date_UTC', 'BY_GSM', electro_property+'_Eastward1', 
                                  electro_property+'_Westward1', 
                                  'Peak_Value_Westward1', 'Peak_Value_Eastward1'])
        df= df[df.Date_UTC.dt.month.isin(months)]
        df=df.replace(9.999e3, np.nan).dropna()
        df_By_pos= df[(df.BY_GSM>=0 )&(df.Peak_Value_Eastward1>0.05)]
        df_By_neg= df[(df.BY_GSM<=0 )&(df.Peak_Value_Eastward1>0.05)]
        cut, w= pd.qcut(df_By_pos.BY_GSM, 20, retbins=True)
        h1= df_By_pos.groupby([cut])
        cut, w= pd.qcut(df_By_neg.BY_GSM, 20, retbins=True)
        h2= df_By_neg.groupby([cut])
        df_By_pos= df[(df.BY_GSM>=0) &(df.Peak_Value_Westward1<-0.1)]
        df_By_neg= df[(df.BY_GSM<=0) &(df.Peak_Value_Westward1<-0.1)]
        cut, w= pd.qcut(df_By_pos.BY_GSM, 20, retbins=True)
        h3= df_By_pos.groupby([cut])
        cut, w= pd.qcut(df_By_neg.BY_GSM, 20, retbins=True)
        h4= df_By_neg.groupby([cut])
        #Mean
        tmp_df['Mean_Westward_pos']= h3.mean()[electro_property+'_Westward1'].values
        tmp_df['Mean_Eastward_pos']= h1.mean()[electro_property+'_Eastward1'].values
        tmp_df['Mean_Westward_neg']= h4.mean()[electro_property+'_Westward1'].values
        tmp_df['Mean_Eastward_neg']= h2.mean()[electro_property+'_Eastward1'].values
        #BY median
        tmp_df['BY_pos_East']= h1.median().BY_GSM.values
        tmp_df['BY_neg_East']= h2.median().BY_GSM.values
        tmp_df['BY_pos_West']= h3.median().BY_GSM.values
        tmp_df['BY_neg_West']= h4.median().BY_GSM.values
        #Errors
        tmp_df['Error_Westward_pos']= h3.std()[electro_property+'_Westward1'].values/np.sqrt(h3.count()[electro_property+'_Westward1']).values
        tmp_df['Error_Eastward_pos']= h1.std()[electro_property+'_Eastward1'].values/np.sqrt(h1.count()[electro_property+'_Eastward1']).values
        tmp_df['Error_Westward_neg']= h4.std()[electro_property+'_Westward1'].values/np.sqrt(h4.count()[electro_property+'_Westward1']).values
        tmp_df['Error_Eastward_neg']= h2.std()[electro_property+'_Eastward1'].values/np.sqrt(h2.count()[electro_property+'_Eastward1']).values
        tmp_df.to_hdf(outpath,key=key,mode='a',append=True,format='t', data_columns=True)
    #Seasonal MLTs
    for mlt in progressbar(np.arange(0, 24), max_value=24):
        mlt1= mlt-1
        mlt2= mlt+1
        if mlt1<0:
            where= "& (MLT>=23 | MLT<=1)"
        else:
            where= "& MLT>=mlt1 & MLT<=mlt2"
        for months, season in zip(np.array([[5,6,7], [11, 12, 1]]), np.array(['summer', 'winter'])):
            key= electro_property+'_'+season+f'_{mlt}'
            if np.any(keys=='/'+key):
                continue
            tmp_df= pd.DataFrame()
            df= store.select(key='main', where="BY_GSM!=9.999e3 & Peak_Value_Westward1!=9.999e3 & Poleward_Boundary_Westward1<UL & Equatorward_Boundary_Westward1>LL"+where,
                             columns=['Date_UTC', 'BY_GSM', electro_property+'_Eastward1', 
                                      electro_property+'_Westward1', 
                                      'Peak_Value_Westward1', 'Peak_Value_Eastward1'])
            df= df[df.Date_UTC.dt.month.isin(months)]
            df=df.replace(9.999e3, np.nan).dropna()
            df_By_pos= df[(df.BY_GSM>=0 )&(df.Peak_Value_Eastward1>0.05)]
            df_By_neg= df[(df.BY_GSM<=0 )&(df.Peak_Value_Eastward1>0.05)]
            cut, w= pd.qcut(df_By_pos.BY_GSM, 20, retbins=True)
            h1= df_By_pos.groupby([cut])
            cut, w= pd.qcut(df_By_neg.BY_GSM, 20, retbins=True)
            h2= df_By_neg.groupby([cut])
            df_By_pos= df[(df.BY_GSM>=0) &(df.Peak_Value_Westward1<-0.1)]
            df_By_neg= df[(df.BY_GSM<=0) &(df.Peak_Value_Westward1<-0.1)]
            cut, w= pd.qcut(df_By_pos.BY_GSM, 20, retbins=True)
            h3= df_By_pos.groupby([cut])
            cut, w= pd.qcut(df_By_neg.BY_GSM, 20, retbins=True)
            h4= df_By_neg.groupby([cut])
            #Mean
            tmp_df['Mean_Westward_pos']= h3.mean()[electro_property+'_Westward1'].values
            tmp_df['Mean_Eastward_pos']= h1.mean()[electro_property+'_Eastward1'].values
            tmp_df['Mean_Westward_neg']= h4.mean()[electro_property+'_Westward1'].values
            tmp_df['Mean_Eastward_neg']= h2.mean()[electro_property+'_Eastward1'].values
            #BY median
            tmp_df['BY_pos_East']= h1.median().BY_GSM.values
            tmp_df['BY_neg_East']= h2.median().BY_GSM.values
            tmp_df['BY_pos_West']= h3.median().BY_GSM.values
            tmp_df['BY_neg_West']= h4.median().BY_GSM.values
            #Errors
            tmp_df['Error_Westward_pos']= h3.std()[electro_property+'_Westward1'].values/np.sqrt(h3.count()[electro_property+'_Westward1']).values
            tmp_df['Error_Eastward_pos']= h1.std()[electro_property+'_Eastward1'].values/np.sqrt(h1.count()[electro_property+'_Eastward1']).values
            tmp_df['Error_Westward_neg']= h4.std()[electro_property+'_Westward1'].values/np.sqrt(h4.count()[electro_property+'_Westward1']).values
            tmp_df['Error_Eastward_neg']= h2.std()[electro_property+'_Eastward1'].values/np.sqrt(h2.count()[electro_property+'_Eastward1']).values
            tmp_df.to_hdf(outpath,key=key,mode='a',append=True,format='t', data_columns=True)
    #All Seasons MLTs
    for mlt in progressbar(np.arange(0, 24), max_value=24):
        mlt1= mlt-1
        mlt2= mlt+1
        if mlt1<0:
            where= "& (MLT>=23 | MLT<=1)"
        else:
            where= "& MLT>=mlt1 & MLT<=mlt2"
        key= electro_property+f'_{mlt}'
        if np.any(keys=='/'+key):
            continue
        tmp_df= pd.DataFrame()
        df= store.select(key='main', where="BY_GSM!=9.999e3 & Peak_Value_Westward1!=9.999e3 & Poleward_Boundary_Westward1<UL & Equatorward_Boundary_Westward1>LL"+where,
                         columns=['Date_UTC', 'BY_GSM', electro_property+'_Eastward1', 
                                  electro_property+'_Westward1', 
                                  'Peak_Value_Westward1', 'Peak_Value_Eastward1'])
        df=df.replace(9.999e3, np.nan).dropna()
        df_By_pos= df[(df.BY_GSM>=0 )&(df.Peak_Value_Eastward1>0.05)]
        df_By_neg= df[(df.BY_GSM<=0 )&(df.Peak_Value_Eastward1>0.05)]
        cut, w= pd.qcut(df_By_pos.BY_GSM, 20, retbins=True)
        h1= df_By_pos.groupby([cut])
        cut, w= pd.qcut(df_By_neg.BY_GSM, 20, retbins=True)
        h2= df_By_neg.groupby([cut])
        df_By_pos= df[(df.BY_GSM>=0) &(df.Peak_Value_Westward1<-0.1)]
        df_By_neg= df[(df.BY_GSM<=0) &(df.Peak_Value_Westward1<-0.1)]
        cut, w= pd.qcut(df_By_pos.BY_GSM, 20, retbins=True)
        h3= df_By_pos.groupby([cut])
        cut, w= pd.qcut(df_By_neg.BY_GSM, 20, retbins=True)
        h4= df_By_neg.groupby([cut])
        #Mean
        tmp_df['Mean_Westward_pos']= h3.mean()[electro_property+'_Westward1'].values
        tmp_df['Mean_Eastward_pos']= h1.mean()[electro_property+'_Eastward1'].values
        tmp_df['Mean_Westward_neg']= h4.mean()[electro_property+'_Westward1'].values
        tmp_df['Mean_Eastward_neg']= h2.mean()[electro_property+'_Eastward1'].values
        #BY median
        tmp_df['BY_pos_East']= h1.median().BY_GSM.values
        tmp_df['BY_neg_East']= h2.median().BY_GSM.values
        tmp_df['BY_pos_West']= h3.median().BY_GSM.values
        tmp_df['BY_neg_West']= h4.median().BY_GSM.values
        #Errors
        tmp_df['Error_Westward_pos']= h3.std()[electro_property+'_Westward1'].values/np.sqrt(h3.count()[electro_property+'_Westward1']).values
        tmp_df['Error_Eastward_pos']= h1.std()[electro_property+'_Eastward1'].values/np.sqrt(h1.count()[electro_property+'_Eastward1']).values
        tmp_df['Error_Westward_neg']= h4.std()[electro_property+'_Westward1'].values/np.sqrt(h4.count()[electro_property+'_Westward1']).values
        tmp_df['Error_Eastward_neg']= h2.std()[electro_property+'_Eastward1'].values/np.sqrt(h2.count()[electro_property+'_Eastward1']).values
        tmp_df.to_hdf(outpath,key=key,mode='a',append=True,format='t', data_columns=True)
    #ALl Seasons All MLTS
    key= electro_property
    if not np.any(keys=='/'+key):
        tmp_df= pd.DataFrame()
        df= store.select(key='main', where="BY_GSM!=9.999e3 & Peak_Value_Westward1!=9.999e3 & Poleward_Boundary_Westward1<UL & Equatorward_Boundary_Westward1>LL",
                         columns=['Date_UTC', 'BY_GSM', electro_property+'_Eastward1', 
                                  electro_property+'_Westward1', 
                                  'Peak_Value_Westward1', 'Peak_Value_Eastward1'])
        df=df.replace(9.999e3, np.nan).dropna()
        df_By_pos= df[(df.BY_GSM>=0 )&(df.Peak_Value_Eastward1>0.05)]
        df_By_neg= df[(df.BY_GSM<=0 )&(df.Peak_Value_Eastward1>0.05)]
        cut, w= pd.qcut(df_By_pos.BY_GSM, 20, retbins=True)
        h1= df_By_pos.groupby([cut])
        cut, w= pd.qcut(df_By_neg.BY_GSM, 20, retbins=True)
        h2= df_By_neg.groupby([cut])
        df_By_pos= df[(df.BY_GSM>=0) &(df.Peak_Value_Westward1<-0.1)]
        df_By_neg= df[(df.BY_GSM<=0) &(df.Peak_Value_Westward1<-0.1)]
        cut, w= pd.qcut(df_By_pos.BY_GSM, 20, retbins=True)
        h3= df_By_pos.groupby([cut])
        cut, w= pd.qcut(df_By_neg.BY_GSM, 20, retbins=True)
        h4= df_By_neg.groupby([cut])
        #Mean
        tmp_df['Mean_Westward_pos']= h3.mean()[electro_property+'_Westward1'].values
        tmp_df['Mean_Eastward_pos']= h1.mean()[electro_property+'_Eastward1'].values
        tmp_df['Mean_Westward_neg']= h4.mean()[electro_property+'_Westward1'].values
        tmp_df['Mean_Eastward_neg']= h2.mean()[electro_property+'_Eastward1'].values
        #BY median
        tmp_df['BY_pos_East']= h1.median().BY_GSM.values
        tmp_df['BY_neg_East']= h2.median().BY_GSM.values
        tmp_df['BY_pos_West']= h3.median().BY_GSM.values
        tmp_df['BY_neg_West']= h4.median().BY_GSM.values
        #Errors
        tmp_df['Error_Westward_pos']= h3.std()[electro_property+'_Westward1'].values/np.sqrt(h3.count()[electro_property+'_Westward1']).values
        tmp_df['Error_Eastward_pos']= h1.std()[electro_property+'_Eastward1'].values/np.sqrt(h1.count()[electro_property+'_Eastward1']).values
        tmp_df['Error_Westward_neg']= h4.std()[electro_property+'_Westward1'].values/np.sqrt(h4.count()[electro_property+'_Westward1']).values
        tmp_df['Error_Eastward_neg']= h2.std()[electro_property+'_Eastward1'].values/np.sqrt(h2.count()[electro_property+'_Eastward1']).values
        tmp_df.to_hdf(outpath,key=key,mode='a',append=True,format='t', data_columns=True)       
    store.close()
    print('>'*5 +f'{electro_property} Complete'+ '<'*5)
    
def profile_create(inpath, outpath):
    import warnings
    warnings.filterwarnings("ignore")
    print('\n','Profile', '\n')
    mlat= np.linspace(49, 81, 50)
    columns= [str(title)+'_Current_East' for title in mlat]
    # columns+=[str(title)+'_Current_North' for title in mlat]
    # columns+=[str(title)+'_Br' for title in mlat]
    store=pd.HDFStore(inpath, mode='r')
    keys=[]
    if os.path.isfile(outpath):
        out_store= pd.HDFStore(outpath, mode='r')
        keys= out_store.keys()
        out_store.close()
    keys=np.array(keys, dtype='object')
    #All Seasons All IMF All MLT All IMF
    key= 'Profile'
    if not np.any(keys=='/'+key):
        df= store.select(key='main', columns= columns+['BY_GSM'], where= "BY_GSM!=9.999e3")
        df=df.replace(9.999e3, np.nan).dropna()
        index_e= np.any(df[columns]>0.05, axis=1)
        index_w= np.any(df[columns]<-0.1, axis=1)
        index_pos= df.BY_GSM>=0
        index_neg= df.BY_GSM<=0
        #East
        cut, w= pd.qcut(abs(df.BY_GSM[index_e]), 3, retbins=True)
        cut= pd.cut(df.BY_GSM[index_e&index_pos], w)
        h= df[index_e&index_pos].groupby(cut)
        cut= pd.cut(abs(df.BY_GSM[index_e&index_neg]), w)
        h2= df[index_e&index_neg].groupby(cut)
        #BY median
        BY_GSM_pos= h.median().BY_GSM.values
        BY_GSM_neg= h2.median().BY_GSM.values
        #Mean
        Currents_pos= h.mean()[columns].values
        Currents_neg= h2.mean()[columns].values
        #Error
        error_pos=h.std()[columns].values/np.sqrt(h.count()[columns]).values
        error_neg=h2.std()[columns].values/np.sqrt(h.count()[columns]).values
        #West
        cut, w= pd.qcut(abs(df.BY_GSM[index_w]), 3, retbins=True)
        cut= pd.cut(df.BY_GSM[index_w&index_pos], w)
        h= df[index_w&index_pos].groupby(cut)
        cut= pd.cut(abs(df.BY_GSM[index_w&index_neg]), w)
        h2= df[index_w&index_neg].groupby(cut)
        #BY median
        BY_GSM_pos= np.append(BY_GSM_pos, h.median().BY_GSM.values)
        BY_GSM_neg= np.append(BY_GSM_neg, h2.median().BY_GSM.values)
        #Mean
        Currents_pos= np.append(Currents_pos, h.mean()[columns].values, axis=0)
        Currents_neg= np.append(Currents_neg, h2.mean()[columns].values, axis=0)
        #Error
        error_pos=np.append(error_pos, h.std()[columns].values/np.sqrt(h.count()[columns]).values, axis=0)
        error_neg=np.append(error_neg, h2.std()[columns].values/np.sqrt(h.count()[columns]).values, axis=0)
        #Output
        d={}
        d.update({columns[i]:np.append(Currents_pos, Currents_neg, axis=0)[:,i] for i in range(50)})
        d.update({'BY_GSM': np.append(BY_GSM_pos, BY_GSM_neg)})
        d.update({f'Error{i}':np.append(error_pos, error_neg, axis=0)[:,i] for i in range(50)})
        d.update({'Jet_BY':['East_pos']*3 + ['West_pos']*3 + ['East_neg']*3 +['West_neg']*3})
        new_df= pd.DataFrame(d)
        # new_df[columns]= np.append(Currents_pos, Currents_neg, axis=0).T
        # new_df['BY_GSM']= np.append(BY_GSM_pos, BY_GSM_neg, axis=0).T
        # new_df['Error']= np.append(error_pos, error_neg, axis=0).T
        # new_df['Jet_BY']= ['East_pos']*20 + ['West_pos']*20 + ['East_neg']*20 +['West_pos']*20
        new_df.to_hdf(outpath,key=key,mode='a',append=True,format='t', data_columns=True)
    #Seasonal All MLTs BY Seperated
    for months, season in progressbar(zip(np.array([[5,6,7], [11, 12, 1]]), np.array(['summer', 'winter'])), max_value=2):
        key= 'Profile'+'_'+season
        if np.any(keys=='/'+key):
            continue
        df= store.select(key='main', columns= columns+['BY_GSM', 'Date_UTC'], where= "BY_GSM!=9.999e3")
        df= df[df.Date_UTC.dt.month.isin(months)]
        df=df.replace(9.999e3, np.nan).dropna()
        index_e= np.any(df[columns]>0.05, axis=1)
        index_w= np.any(df[columns]<-0.1, axis=1)
        index_pos= df.BY_GSM>=0
        index_neg= df.BY_GSM<=0
        #East
        cut, w= pd.qcut(abs(df.BY_GSM[index_e]), 3, retbins=True)
        cut= pd.cut(df.BY_GSM[index_e&index_pos], w)
        h= df[index_e&index_pos].groupby(cut)
        cut= pd.cut(abs(df.BY_GSM[index_e&index_neg]), w)
        h2= df[index_e&index_neg].groupby(cut)
        #BY median
        BY_GSM_pos= h.median().BY_GSM.values
        BY_GSM_neg= h2.median().BY_GSM.values
        #Mean
        Currents_pos= h.mean()[columns].values
        Currents_neg= h2.mean()[columns].values
        #Error
        error_pos=h.std()[columns].values/np.sqrt(h.count()[columns]).values
        error_neg=h2.std()[columns].values/np.sqrt(h.count()[columns]).values
        #West
        cut, w= pd.qcut(abs(df.BY_GSM[index_w]), 3, retbins=True)
        cut= pd.cut(df.BY_GSM[index_w&index_pos], w)
        h= df[index_w&index_pos].groupby(cut)
        cut= pd.cut(abs(df.BY_GSM[index_w&index_neg]), w)
        h2= df[index_w&index_neg].groupby(cut)
        #BY median
        BY_GSM_pos= np.append(BY_GSM_pos, h.median().BY_GSM.values)
        BY_GSM_neg= np.append(BY_GSM_neg, h2.median().BY_GSM.values)
        #Mean
        Currents_pos= np.append(Currents_pos, h.mean()[columns].values, axis=0)
        Currents_neg= np.append(Currents_neg, h2.mean()[columns].values, axis=0)
        #Error
        error_pos=np.append(error_pos, h.std()[columns].values/np.sqrt(h.count()[columns]).values, axis=0)
        error_neg=np.append(error_neg, h2.std()[columns].values/np.sqrt(h.count()[columns]).values, axis=0)
        #Output
        d={}
        d.update({columns[i]:np.append(Currents_pos, Currents_neg, axis=0)[:,i] for i in range(50)})
        d.update({'BY_GSM': np.append(BY_GSM_pos, BY_GSM_neg)})
        d.update({f'Error{i}':np.append(error_pos, error_neg, axis=0)[:,i] for i in range(50)})
        d.update({'Jet_BY':['East_pos']*3 + ['West_pos']*3 + ['East_neg']*3 +['West_neg']*3})
        new_df= pd.DataFrame(d)
        # new_df[columns]= np.append(Currents_pos, Currents_neg, axis=0).T
        # new_df['BY_GSM']= np.append(BY_GSM_pos, BY_GSM_neg, axis=0).T
        # new_df['Error']= np.append(error_pos, error_neg, axis=0).T
        # new_df['Jet_BY']= ['East_pos']*20 + ['West_pos']*20 + ['East_neg']*20 +['West_pos']*20
        new_df.to_hdf(outpath,key=key,mode='a',append=True,format='t', data_columns=True)
    #All Seasons MLTs BY Seperated
    for mlt in progressbar(np.arange(0,24), max_value=24):
        mlt1=mlt-1
        mlt2=mlt+1
        if mlt1<0:
            where= "& (MLT>=23 | MLT<=1)"
        else:
            where= "& MLT>=mlt1 & MLT<=mlt2"
        key= f'Profile_{mlt}'
        if np.any(keys=='/'+key):
            continue
        df= store.select(key='main', columns= columns+['BY_GSM'], where= "BY_GSM!=9.999e3"+where)
        df=df.replace(9.999e3, np.nan).dropna()
        index_e= np.any(df[columns]>0.05, axis=1)
        index_w= np.any(df[columns]<-0.1, axis=1)
        index_pos= df.BY_GSM>=0
        index_neg= df.BY_GSM<=0
        #East
        cut, w= pd.qcut(abs(df.BY_GSM[index_e]), 3, retbins=True)
        cut= pd.cut(df.BY_GSM[index_e&index_pos], w)
        h= df[index_e&index_pos].groupby(cut)
        cut= pd.cut(abs(df.BY_GSM[index_e&index_neg]), w)
        h2= df[index_e&index_neg].groupby(cut)
        #BY median
        BY_GSM_pos= h.median().BY_GSM.values
        BY_GSM_neg= h2.median().BY_GSM.values
        #Mean
        Currents_pos= h.mean()[columns].values
        Currents_neg= h2.mean()[columns].values
        #Error
        error_pos=h.std()[columns].values/np.sqrt(h.count()[columns]).values
        error_neg=h2.std()[columns].values/np.sqrt(h.count()[columns]).values
        #West
        cut, w= pd.qcut(abs(df.BY_GSM[index_w]), 3, retbins=True)
        cut= pd.cut(df.BY_GSM[index_w&index_pos], w)
        h= df[index_w&index_pos].groupby(cut)
        cut= pd.cut(abs(df.BY_GSM[index_w&index_neg]), w)
        h2= df[index_w&index_neg].groupby(cut)
        #BY median
        BY_GSM_pos= np.append(BY_GSM_pos, h.median().BY_GSM.values)
        BY_GSM_neg= np.append(BY_GSM_neg, h2.median().BY_GSM.values)
        #Mean
        Currents_pos= np.append(Currents_pos, h.mean()[columns].values, axis=0)
        Currents_neg= np.append(Currents_neg, h2.mean()[columns].values, axis=0)
        #Error
        error_pos=np.append(error_pos, h.std()[columns].values/np.sqrt(h.count()[columns]).values, axis=0)
        error_neg=np.append(error_neg, h2.std()[columns].values/np.sqrt(h.count()[columns]).values, axis=0)
        #Output
        d={}
        d.update({columns[i]:np.append(Currents_pos, Currents_neg, axis=0)[:,i] for i in range(50)})
        d.update({'BY_GSM': np.append(BY_GSM_pos, BY_GSM_neg)})
        d.update({f'Error{i}':np.append(error_pos, error_neg, axis=0)[:,i] for i in range(50)})
        d.update({'Jet_BY':['East_pos']*3 + ['West_pos']*3 + ['East_neg']*3 +['West_neg']*3})
        new_df= pd.DataFrame(d)
        # new_df[columns]= np.append(Currents_pos, Currents_neg, axis=0).T
        # new_df['BY_GSM']= np.append(BY_GSM_pos, BY_GSM_neg, axis=0).T
        # new_df['Error']= np.append(error_pos, error_neg, axis=0).T
        # new_df['Jet_BY']= ['East_pos']*20 + ['West_pos']*20 + ['East_neg']*20 +['West_pos']*20
        new_df.to_hdf(outpath,key=key,mode='a',append=True,format='t', data_columns=True)
        
    # Seasonal MLT BY seperated
    for mlt in progressbar(np.arange(0, 24), max_value=24):
        mlt1= mlt-1
        mlt2= mlt+1
        if mlt1<0:
            where= "& (MLT>=23 | MLT<=1)"
        else:
            where= "& MLT>=mlt1 & MLT<=mlt2"
        for months, season in zip(np.array([[5,6,7], [11, 12, 1]]), np.array(['summer', 'winter'])):
            key= 'Profile'+'_'+season+f'_{mlt}'
            if np.any(keys=='/'+key):
                continue
            df= store.select(key='main', columns= columns+['BY_GSM', 'Date_UTC'], where= "BY_GSM!=9.999e3"+where)
            df= df[df.Date_UTC.dt.month.isin(months)]
            df=df.replace(9.999e3, np.nan).dropna()
            index_e= np.any(df[columns]>0.05, axis=1)
            index_w= np.any(df[columns]<-0.1, axis=1)
            index_pos= df.BY_GSM>=0
            index_neg= df.BY_GSM<=0
            #East
            cut, w= pd.qcut(abs(df.BY_GSM[index_e]), 3, retbins=True)
            cut= pd.cut(df.BY_GSM[index_e&index_pos], w)
            h= df[index_e&index_pos].groupby(cut)
            cut= pd.cut(abs(df.BY_GSM[index_e&index_neg]), w)
            h2= df[index_e&index_neg].groupby(cut)
            #BY median
            BY_GSM_pos= h.median().BY_GSM.values
            BY_GSM_neg= h2.median().BY_GSM.values
            #Mean
            Currents_pos= h.mean()[columns].values
            Currents_neg= h2.mean()[columns].values
            #Error
            error_pos=h.std()[columns].values/np.sqrt(h.count()[columns]).values
            error_neg=h2.std()[columns].values/np.sqrt(h.count()[columns]).values
            #West
            cut, w= pd.qcut(abs(df.BY_GSM[index_w]), 3, retbins=True)
            cut= pd.cut(df.BY_GSM[index_w&index_pos], w)
            h= df[index_w&index_pos].groupby(cut)
            cut= pd.cut(abs(df.BY_GSM[index_w&index_neg]), w)
            h2= df[index_w&index_neg].groupby(cut)
            #BY median
            BY_GSM_pos= np.append(BY_GSM_pos, h.median().BY_GSM.values)
            BY_GSM_neg= np.append(BY_GSM_neg, h2.median().BY_GSM.values)
            #Mean
            Currents_pos= np.append(Currents_pos, h.mean()[columns].values, axis=0)
            Currents_neg= np.append(Currents_neg, h2.mean()[columns].values, axis=0)
            #Error
            error_pos=np.append(error_pos, h.std()[columns].values/np.sqrt(h.count()[columns]).values, axis=0)
            error_neg=np.append(error_neg, h2.std()[columns].values/np.sqrt(h.count()[columns]).values, axis=0)
            #Output
            d={}
            d.update({columns[i]:np.append(Currents_pos, Currents_neg, axis=0)[:,i] for i in range(50)})
            d.update({'BY_GSM': np.append(BY_GSM_pos, BY_GSM_neg)})
            d.update({f'Error{i}':np.append(error_pos, error_neg, axis=0)[:,i] for i in range(50)})
            d.update({'Jet_BY':['East_pos']*3 + ['West_pos']*3 + ['East_neg']*3 +['West_neg']*3})
            new_df= pd.DataFrame(d)
            # new_df[columns]= np.append(Currents_pos, Currents_neg, axis=0).T
            # new_df['BY_GSM']= np.append(BY_GSM_pos, BY_GSM_neg, axis=0).T
            # new_df['Error']= np.append(error_pos, error_neg, axis=0).T
            # new_df['Jet_BY']= ['East_pos']*20 + ['West_pos']*20 + ['East_neg']*20 +['West_pos']*20
            new_df.to_hdf(outpath,key=key,mode='a',append=True,format='t', data_columns=True)
    #All Seasons Clock Angle All MLTs
    for Clock_Angle in progressbar(np.arange(-180, 190, 45), max_value=9):
        key= 'Profile'+f'_Clock:{Clock_Angle}'
        if np.any(keys=='/'+key):
            continue
        if abs(Clock_Angle)== 180:
            where= "& (Clock_Angle>= 157.5 | Clock_Angle<=-157.5)"
        else:
            clock1= Clock_Angle-22.5
            clock2= Clock_Angle+22.5
            where= "& Clock_Angle>=clock1 & Clock_Angle<=clock2"
        df= store.select(key='main', columns= columns+['BY_GSM', 'Date_UTC'], where= "BY_GSM!=9.999e3"+where)
        df=df.replace(9.999e3, np.nan).dropna()
        index_e= np.any(df[columns]>0.05, axis=1)
        index_w= np.any(df[columns]<-0.1, axis=1)
        #Mean
        Currents= df[index_e].mean()[columns].values
        #Error
        error=df[index_e].std()[columns].values/np.sqrt(df[index_e].count()[columns]).values
        #West
        #Mean
        Currents= np.vstack([Currents, df[index_w].mean()[columns].values])
        #Error
        error=np.vstack([error, df[index_w].std()[columns].values/np.sqrt(df[index_e].count()[columns]).values])
        error=error.astype('float64')
        #Output
        d={}
        d.update({columns[i]:Currents[:,i] for i in range(50)})
        d.update({f'Error{i}':error[:,i] for i in range(50)})
        d.update({'Jet':['East'] + ['West']})
        new_df= pd.DataFrame(d)
        # new_df[columns]= np.append(Currents_pos, Currents_neg, axis=0).T
        # new_df['BY_GSM']= np.append(BY_GSM_pos, BY_GSM_neg, axis=0).T
        # new_df['Error']= np.append(error_pos, error_neg, axis=0).T
        # new_df['Jet_BY']= ['East_pos']*20 + ['West_pos']*20 + ['East_neg']*20 +['West_pos']*20
        new_df.to_hdf(outpath,key=key,mode='a',append=True,format='t', data_columns=True)
    #All Seasons Clock Angle MLTs
    for mlt in progressbar(np.arange(0, 24), max_value=24):
        mlt1= mlt-1
        mlt2= mlt+1
        if mlt1<0:
            where= "& (MLT>=23 | MLT<=1)"
        else:
            where= "& MLT>=mlt1 & MLT<=mlt2"
        for Clock_Angle in progressbar(np.arange(-180, 190, 45), max_value=9):
            key= 'Profile'+f'_{mlt}_Clock:{Clock_Angle}'
            if np.any(keys=='/'+key):
                continue
            if abs(Clock_Angle)== 180:
                where2= "& (Clock_Angle>= 157.5 | Clock_Angle<=-157.5)"
            else:
                clock1= Clock_Angle-22.5
                clock2= Clock_Angle+22.5
                where2= "& Clock_Angle>=clock1 & Clock_Angle<=clock2"
            df= store.select(key='main', columns= columns+['BY_GSM', 'Date_UTC'], where= "BY_GSM!=9.999e3"+where+where2)
            df=df.replace(9.999e3, np.nan).dropna()
            index_e= np.any(df[columns]>0.05, axis=1)
            index_w= np.any(df[columns]<-0.1, axis=1)
            #Mean
            Currents= df[index_e].mean()[columns].values
            #Error
            error=df[index_e].std()[columns].values/np.sqrt(df[index_e].count()[columns]).values
            #West
            #Mean
            Currents= np.vstack([Currents, df[index_w].mean()[columns].values])
            #Error
            error=np.vstack([error, df[index_w].std()[columns].values/np.sqrt(df[index_e].count()[columns]).values])
            error=error.astype('float64')
            #Output
            d={}
            d.update({columns[i]:Currents[:,i] for i in range(50)})
            d.update({f'Error{i}':error[:,i] for i in range(50)})
            d.update({'Jet':['East'] + ['West']})
            new_df= pd.DataFrame(d)
            # new_df[columns]= np.append(Currents_pos, Currents_neg, axis=0).T
            # new_df['BY_GSM']= np.append(BY_GSM_pos, BY_GSM_neg, axis=0).T
            # new_df['Error']= np.append(error_pos, error_neg, axis=0).T
            # new_df['Jet_BY']= ['East_pos']*20 + ['West_pos']*20 + ['East_neg']*20 +['West_pos']*20
            new_df.to_hdf(outpath,key=key,mode='a',append=True,format='t', data_columns=True)
    #Seasonal Clock Angle All MLTS
    for Clock_Angle in progressbar(np.arange(-180, 190, 45), max_value=9):
        if abs(Clock_Angle)== 180:
            where= "& (Clock_Angle>= 157.5 | Clock_Angle<=-157.5)"
        else:
            clock1= Clock_Angle-22.5
            clock2= Clock_Angle+22.5
            where= "& Clock_Angle>=clock1 & Clock_Angle<=clock2"
        for months, season in zip(np.array([[5,6,7], [11, 12, 1]]), np.array(['summer', 'winter'])):
            key= 'Profile'+f'_{season}_Clock:{Clock_Angle}'
            if np.any(keys=='/'+key):
                continue            
            df= store.select(key='main', columns= columns+['BY_GSM', 'Date_UTC'], where= "BY_GSM!=9.999e3"+where)
            df= df[df.Date_UTC.dt.month.isin(months)]
            df=df.replace(9.999e3, np.nan).dropna()
            index_e= np.any(df[columns]>0.05, axis=1)
            index_w= np.any(df[columns]<-0.1, axis=1)
            #Mean
            Currents= df[index_e].mean()[columns].values
            #Error
            error=df[index_e].std()[columns].values/np.sqrt(df[index_e].count()[columns]).values
            #West
            #Mean
            Currents= np.vstack([Currents, df[index_w].mean()[columns].values])
            #Error
            error=np.vstack([error, df[index_w].std()[columns].values/np.sqrt(df[index_e].count()[columns]).values])
            error=error.astype('float64')
            #Output
            d={}
            d.update({columns[i]:Currents[:,i] for i in range(50)})
            d.update({f'Error{i}':error[:,i] for i in range(50)})
            d.update({'Jet':['East'] + ['West']})
            new_df= pd.DataFrame(d)
            # new_df[columns]= np.append(Currents_pos, Currents_neg, axis=0).T
            # new_df['BY_GSM']= np.append(BY_GSM_pos, BY_GSM_neg, axis=0).T
            # new_df['Error']= np.append(error_pos, error_neg, axis=0).T
            # new_df['Jet_BY']= ['East_pos']*20 + ['West_pos']*20 + ['East_neg']*20 +['West_pos']*20
            new_df.to_hdf(outpath,key=key,mode='a',append=True,format='t', data_columns=True)
    #Seasonal Clock Angle MLTs
    for mlt in progressbar(np.arange(0, 24), max_value=24):
        mlt1= mlt-1
        mlt2= mlt+1
        if mlt1<0:
            where= "& (MLT>=23 | MLT<=1)"
        else:
            where= "& MLT>=mlt1 & MLT<=mlt2"
        for Clock_Angle in progressbar(np.arange(-180, 190, 45), max_value=9):
            if abs(Clock_Angle)== 180:
                where2= "& (Clock_Angle>= 157.5 | Clock_Angle<=-157.5)"
            else:
                clock1= Clock_Angle-22.5
                clock2= Clock_Angle+22.5
                where2= "& Clock_Angle>=clock1 & Clock_Angle<=clock2"
            for months, season in zip(np.array([[5,6,7], [11, 12, 1]]), np.array(['summer', 'winter'])):
                key= 'Profile'+f'_{season}_{mlt}_Clock:{Clock_Angle}'
                if np.any(keys=='/'+key):
                    continue            
                df= store.select(key='main', columns= columns+['BY_GSM', 'Date_UTC'], where= "BY_GSM!=9.999e3"+where+where2)
                df= df[df.Date_UTC.dt.month.isin(months)]
                df=df.replace(9.999e3, np.nan).dropna()
                index_e= np.any(df[columns]>0.05, axis=1)
                index_w= np.any(df[columns]<-0.1, axis=1)
                #Mean
                Currents= df[index_e].mean()[columns].values
                #Error
                error=df[index_e].std()[columns].values/np.sqrt(df[index_e].count()[columns]).values
                #West
                #Mean
                Currents= np.vstack([Currents, df[index_w].mean()[columns].values])
                #Error
                error=np.vstack([error, df[index_w].std()[columns].values/np.sqrt(df[index_e].count()[columns]).values])
                error=error.astype('float64')
                #Output
                d={}
                d.update({columns[i]:Currents[:,i] for i in range(50)})
                d.update({f'Error{i}':error[:,i] for i in range(50)})
                d.update({'Jet':['East'] + ['West']})
                new_df= pd.DataFrame(d)
                # new_df[columns]= np.append(Currents_pos, Currents_neg, axis=0).T
                # new_df['BY_GSM']= np.append(BY_GSM_pos, BY_GSM_neg, axis=0).T
                # new_df['Error']= np.append(error_pos, error_neg, axis=0).T
                # new_df['Jet_BY']= ['East_pos']*20 + ['West_pos']*20 + ['East_neg']*20 +['West_pos']*20
                new_df.to_hdf(outpath,key=key,mode='a',append=True,format='t', data_columns=True)
    store.close()
    print('>'*5 +'Profiles Complete'+ '<'*5)
def average_current_create(inpath, outpath):
    import warnings
    warnings.filterwarnings("ignore")
    print('\n','Average Current', '\n')
    mlat= np.linspace(49, 81, 50)
    columns= [str(title)+'_Current_East' for title in mlat]
    columns+=[str(title)+'_Current_North' for title in mlat]
    columns+=[str(title)+'_Br' for title in mlat]
    store=pd.HDFStore(inpath, mode='r')
    keys=[]
    if os.path.isfile(outpath):
        out_store= pd.HDFStore(outpath, mode='r')
        keys= out_store.keys()
        out_store.close()
    keys=np.array(keys, dtype='object')
    #All
    key= 'Average_Current'
    if np.all(keys!='/'+key):
        df= store.select(key='main', columns= columns+['Date_UTC', 'MLT'])
        df=df.replace(9.999e3, np.nan).dropna()
        cut= pd.cut(df.MLT, np.arange(0, 24.5, .5))
        h= df.groupby([cut])
        Current_east= h.mean()[columns[:50]].values
        Current_north= h.mean()[columns[50:100]].values
        Br= h.mean()[columns[100:]].values
        d={f'East_{mlt}':Current_east[i] for i, mlt in enumerate(np.arange(0.25, 24, .5))}
        d.update({f'North_{mlt}':Current_north[i] for i, mlt in enumerate(np.arange(0.25, 24, .5))})
        d.update({f'Br_{mlt}':Br[i] for i, mlt in enumerate(np.arange(0.25, 24, .5))})
        new_df= pd.DataFrame(d)
        new_df.to_hdf(outpath, key=key, mode='a', append=True, format='t', data_columns=True)
    #Seasonal
    for months, season in progressbar(zip(np.array([[5,6,7], [11, 12, 1]]), np.array(['summer', 'winter'])), max_value=2):
        key= 'Average_Current'+'_'+season
        if np.any(keys=='/'+key):
            continue
        df= store.select(key='main', columns= columns+['Date_UTC', 'MLT'])
        df= df[df.Date_UTC.dt.month.isin(months)]
        df=df.replace(9.999e3, np.nan).dropna()
        cut= pd.cut(df.MLT, np.arange(0, 24.5, .5))
        h= df.groupby([cut])
        Current_east= h.mean()[columns[:50]].values
        Current_north= h.mean()[columns[50:100]].values
        Br= h.mean()[columns[100:]].values
        d={f'East_{mlt}':Current_east[i] for i, mlt in enumerate(np.arange(0.25, 24, .5))}
        d.update({f'North_{mlt}':Current_north[i] for i, mlt in enumerate(np.arange(0.25, 24, .5))})
        d.update({f'Br_{mlt}':Br[i] for i, mlt in enumerate(np.arange(0.25, 24, .5))})
        new_df= pd.DataFrame(d)
        new_df.to_hdf(outpath, key=key, mode='a', append=True, format='t', data_columns=True)
    #All Seasons Clock Angle
    for Clock_Angle in progressbar(np.arange(-180, 190, 45), max_value=9):
        key= 'Average_Current'+f'_Clock:{Clock_Angle}'
        if np.any(keys=='/'+key):
            continue
        if abs(Clock_Angle)== 180:
            where= "(Clock_Angle>= 157.5 | Clock_Angle<=-157.5)"
        else:
            clock1= Clock_Angle-22.5
            clock2= Clock_Angle+22.5
            where= "Clock_Angle>=clock1 & Clock_Angle<=clock2"
        df= store.select(key='main', columns= columns+['Date_UTC', 'MLT'], where= where)
        df=df.replace(9.999e3, np.nan).dropna()
        cut= pd.cut(df.MLT, np.arange(0, 24.5, .5))
        h= df.groupby([cut])
        Current_east= h.mean()[columns[:50]].values
        Current_north= h.mean()[columns[50:100]].values
        Br= h.mean()[columns[100:]].values
        d={f'East_{mlt}':Current_east[i] for i, mlt in enumerate(np.arange(0.25, 24, .5))}
        d.update({f'North_{mlt}':Current_north[i] for i, mlt in enumerate(np.arange(0.25, 24, .5))})
        d.update({f'Br_{mlt}':Br[i] for i, mlt in enumerate(np.arange(0.25, 24, .5))})
        new_df= pd.DataFrame(d)
        new_df.to_hdf(outpath, key=key, mode='a', append=True, format='t', data_columns=True)
    #Seasonal Clock Angle
    for months, season in progressbar(zip(np.array([[5,6,7], [11, 12, 1]]), np.array(['summer', 'winter'])), max_value=2):
        for Clock_Angle in progressbar(np.arange(-180, 190, 45), max_value=9):
            key= 'Average_Current'+f'_{season}_Clock:{Clock_Angle}'
            if np.any(keys=='/'+key):
                continue
            if abs(Clock_Angle)== 180:
                where= "(Clock_Angle>= 157.5 | Clock_Angle<=-157.5)"
            else:
                clock1= Clock_Angle-22.5
                clock2= Clock_Angle+22.5
                where= "Clock_Angle>=clock1 & Clock_Angle<=clock2"
            df= store.select(key='main', columns= columns+['Date_UTC', 'MLT'], where= where)
            df= df[df.Date_UTC.dt.month.isin(months)]
            df=df.replace(9.999e3, np.nan).dropna()
            cut= pd.cut(df.MLT, np.arange(0, 24.5, .5))
            h= df.groupby([cut])
            Current_east= h.mean()[columns[:50]].values
            Current_north= h.mean()[columns[50:100]].values
            Br= h.mean()[columns[100:]].values
            d={f'East_{mlt}':Current_east[i] for i, mlt in enumerate(np.arange(0.25, 24, .5))}
            d.update({f'North_{mlt}':Current_north[i] for i, mlt in enumerate(np.arange(0.25, 24, .5))})
            d.update({f'Br_{mlt}':Br[i] for i, mlt in enumerate(np.arange(0.25, 24, .5))})
            new_df= pd.DataFrame(d)
            new_df.to_hdf(outpath, key=key, mode='a', append=True, format='t', data_columns=True)
    store.close()
    print('>'*5 +'Average Currents Complete'+ '<'*5)
def Create(boundaries_path, meridian_path, out_path):
    profile_create(meridian_path, out_path)
    average_current_create(meridian_path, out_path)
if __name__=='__main__':
    # boundaries_path= './SECS_Data/electrojet_boundaries_Clock.hdf5'
    # out_path= './Plotting_Data/SECS_Visualisation.hdf5'
    # meridian_path= './SECS_Data/results_singularity_mod_clock2.hdf5'
    # Create(boundaries_path, meridian_path, out_path)
    boundaries_path='./SECS_Data/electrojet_boundaries_IMF_vaex.hdf5'
    out_path='./Plotting_Data/SECS_Visualisation_Vaex.hdf5'
    boundary_create(boundaries_path, out_path, 'Width')
# boundary_create(boundaries_path, out_path, 'Current')
# boundary_create(boundaries_path, out_path, 'Peak_Value')
# profile_create(meridian_path, out_path)
# average_current_create(meridian_path, out_path)

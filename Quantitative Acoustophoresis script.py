# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 14:24:25 2021

@author: bvadg
"""

import os
import math
import numpy as np
import pandas as pd
import scipy.optimize
import ptitprince as pt
import matplotlib.pyplot as plt
from matplotlib import cm
from nptdms import TdmsFile
from scipy.optimize import curve_fit
from lmfit import Model, Parameter, fit_report, minimize
from scipy.interpolate import RegularGridInterpolator, interp2d, griddata




"""
Parameters and settings
"""
dir_calibr = r"E:\MyFiles\Experimental data\Acoustophoresis\20200923 - PMNs and empty liposomes 1 of 3\AP_4.65um_SAV_coated_PS_beads_10-0p45V_20x_recalibration"
dir_data = r"E:\MyFiles\Experimental data\Acoustophoresis\20200923 - PMNs and empty liposomes 1 of 3\AP_PMNs_HEPES+_10-0p65V_20x\3"
dir_script = r"E:\MyFiles\Scripts\Acoustophoresis"

temperature = 20                    # Celsius degrees; room temperature or the one set using a special stage
objective_magn = 20                 # 10/20/40 depending on which objective is used
diameter_sample = 4.65              # microbead size

is_size_known = False               # True - for monodisperse microbeads, False for other samples
is_temp_dep = False                 # True - use the temperature to calculate precisely the viscosity of the medium
is_density_known = False             # True - the sample density is known, False - otherwise
is_calibration = False              # True - for calibration sample, False - otherwise
is_time_resolved = False
calibration_method = "average"      # "heatmap" or "average"       
density_method = "full"             # "linear", "full" or "forces" fit
sample = "PMN"                      # material of the sample (PS/Silica/PMMA) or its type (MEF/PMN/GUV/GUV gel)
medium = "HEPES+"                     # PBS/water/glucose/HEPES+

corr_bead = 1       # the correction factor that accounts for the different properties of the 5.31 um PS beads from a different supplier


#%%
"""
Creation of the dictionary with the physical parameters
Functions to load the data, parse it and prepare for the analysisfde
"""
    
def set_config_values():
    os.chdir(dir_script)
    file_physical_properties = "Physical properties.xlsx"
    properties_samples = pd.read_excel(file_physical_properties, index_col="Sample", sheet_name="Samples")
    properties_mediums = pd.read_excel(file_physical_properties, index_col="Medium", sheet_name="Mediums")
    
    Phys_properties = {}
    
    density_sample = float(properties_samples.loc[sample]['Density'])
    Phys_properties['density_sample'] = density_sample # kg/m3
    density_sample_std = float(properties_samples.loc[sample]['DensitySTD'])
    Phys_properties['density_sample_std'] = density_sample_std # kg/m3
    density_medium = float(properties_mediums.loc[medium]['Density']) 
    Phys_properties['density_medium'] = density_medium # kg/m3
    Phys_properties['temperature'] = temperature  # oC
    
    Phys_properties['nm_per_px'] = calculate_nm_per_px(objective_magn) # calculated based on the objective in use
    Phys_properties['nm_per_px_def'] = 296 # default value in the LabView software
    Phys_properties['FOV_x_px'] = 1280 # field of view x-size in px
    Phys_properties['FOV_y_px'] = 1024 # field of view y-size in px
    Phys_properties['FOV_x_um'] = Phys_properties['FOV_x_px'] * Phys_properties['nm_per_px'] / 1000
    Phys_properties['FOV_y_um'] = Phys_properties['FOV_y_px'] * Phys_properties['nm_per_px'] / 1000
    
    viscosity_medium_at_20  = float(properties_mediums.loc[medium]['Viscosity'])    
    compr_medium = float(properties_mediums.loc[medium]['Compressibility']) 
    Phys_properties['compr_medium'] = compr_medium # 10^-10 Pa^-1   
    Phys_properties['compr_PS'] = 2.30 # 10^-10 Pa^-1 
    Phys_properties['g'] = 9.813 # m/s^2
    Phys_properties['frequency_resonance'] = 14.5 # MHz
    Phys_properties['velocity_sound'] = 1480 # m/s at 20 degrees Celsius
    
    if is_size_known == True:
        Phys_properties['diameter_sample'] = diameter_sample #um
        Phys_properties['radius_sample'] = diameter_sample / 2 # um
    if is_temp_dep == False:
        Phys_properties['viscosity_medium'] = viscosity_medium_at_20  # Pa*s at 20 degrees Celsius
    else:
        Phys_properties['viscosity_medium'] = viscosity_medium_at_20 * (0.22651 + 1.45718 * np.exp(-temperature / 32.56))    
    if is_calibration == False:
        os.chdir(dir_calibr)
        file_calibration = "Calibration values.xlsx"
        properties_samples = pd.read_excel(file_calibration)
        Phys_properties['radius_calibr'] = properties_samples.iloc[0][0] # um
        Phys_properties['radius_calibr_std'] = properties_samples.iloc[0][1]
        Phys_properties['contrast_factor_calibr'] = properties_samples.iloc[0][2]
        if calibration_method == "average":
            Phys_properties['calibr'] = properties_samples.iloc[0][3] # pN / V^2 response of the calibration sample on that day
            Phys_properties['calibr_std'] = properties_samples.iloc[0][4]
        

    Phys_properties['acoustic_node_position'] = 19 #um
    
    # acoustic node
    # sound velocity
    # calibration size 
    
    return Phys_properties
 
def select_folders():
    if is_calibration == False: # set working directory 
        dir_analysis = dir_data  
        folders = [dir_analysis]
    elif calibration_method == "average":
        dir_analysis = dir_calibr
        folders = [dir_analysis]
    else:
        folders = [f.path for f in os.scandir(dir_calibr) if f.is_dir()]
    return folders
     
def calculate_nm_per_px(magnification):
    if magnification == 10:
        return 1179
    elif magnification == 20:
        return 589
    elif magnification == 40:
        return 294

def load_TVXYZ_data(directory):
    Time = pd.DataFrame()
    Voltage = pd.DataFrame()
    Z = pd.DataFrame()
    XY=pd.DataFrame()
    num=0
    annotation = 'Loading experimental TVXYZ data...'
    print (annotation)
    for filename in os.listdir(directory):
        if filename.endswith('.tdms'):
            RawData = TdmsFile(filename).as_dataframe()    
            RawColumns = RawData.columns
            data = RawData.rename(columns = {RawColumns[0]: 'Time (ms)', 
                                             RawColumns[1]: 'X (um)', 
                                             RawColumns[2]: 'Y (um)', 
                                             RawColumns[3]: 'Z (nm)', 
                                             RawColumns[4]: 'V'})
            if num == 0:
                Time = data['Time (ms)'].copy()
                Voltage = data['V'].copy()
                Z['Z1 (um)'] = data['Z (nm)'].copy()
                XY = data[['X (um)','Y (um)']].copy()
            else:
                Z=pd.concat([Z, data['Z (nm)']], axis = 1)
                XY=pd.concat([XY, data['X (um)'],data['Y (um)']], axis = 1)
            num = num+1
            new_name_x = 'X' + str(num) + ' (um)'
            new_name_y = 'Y' + str(num) + ' (um)'
            new_name_z = 'Z' + str(num) + ' (um)'
            XY.rename(columns = {'X (um)':new_name_x,'Y (um)':new_name_y}, inplace=True)
            Z.rename(columns = {'Z (nm)':new_name_z}, inplace=True)
    XY = XY / 1000 / Phys_properties['nm_per_px_def'] * Phys_properties['nm_per_px']
    XY = XY - XY.iloc[0]
    Z = Z / 1000 
    return Time, Voltage, XY, Z

def load_TVZ_data(directory):
    Time = pd.DataFrame()
    Voltage = pd.DataFrame()
    Z = pd.DataFrame()    
    num = 0
    annotation = 'Loading experimental TVZ data...'
    print (annotation)
    for filename in os.listdir(directory):
        if filename.endswith('.tdms'):
            RawData = TdmsFile(filename).as_dataframe()    
            RawColumns = RawData.columns
            data = RawData.rename(columns = {RawColumns[0]: 'Time (ms)', 
                                             RawColumns[1]: 'X (nm)', 
                                             RawColumns[2]: 'Y (nm)', 
                                             RawColumns[3]: 'Z (nm)', 
                                             RawColumns[4]: 'V'})    
            data.drop(['X (nm)', 'Y (nm)'], axis=1,inplace=True)
            if num==0:
                Time = data['Time (ms)'].copy()
                Voltage = data['V'].copy()
                Z['Z1 (um)'] = data['Z (nm)'].copy()
            else:
                Z=pd.concat([Z, data['Z (nm)']], axis = 1)
            num=num+1
            new_name = 'Z' + str(num) + ' (um)'
            Z.rename(columns={'Z (nm)':new_name}, inplace=True)
    Z = Z / 1000    
    return Time, Voltage, Z

def calculate_xy_positions(XY, directory): 
    ROI=pd.DataFrame()
    starting_XY=pd.DataFrame(columns = ['X0 (um)', 'Y0 (um)'])    
    for filename in os.listdir(directory):
        if filename.endswith('.roi'):
            ROI = pd.read_csv(filename, header=None, sep='\t').transpose()
    ROI.columns = ['X_TopLeft (px)','Y_TopLeft (px)','X_BottomRight (px)','Y_BottomRight (px)']
    starting_XY['X0 (um)']=(ROI['X_TopLeft (px)']+ROI['X_BottomRight (px)']) / 2 * Phys_properties['nm_per_px'] / 1000
    starting_XY['Y0 (um)']=(ROI['Y_TopLeft (px)']+ROI['Y_BottomRight (px)']) / 2 * Phys_properties['nm_per_px'] / 1000
    XY = XY + starting_XY.values.flatten()
    return XY

def find_voltag_spikes(Voltage):
    voltage_map=pd.DataFrame(columns=['Amplitude (V)','Start','End'])
    f = 0
    for count, V in enumerate(Voltage):
        if V > 0 and f == 0:
            voltage_map = voltage_map.append(dict(zip(voltage_map.columns, [V,count-1,0])), ignore_index=True)
            f = 1
        elif V == 0 and f == 1:
            voltage_map.at[voltage_map.index[-1],'End'] = count - 1
            f = 0
    
    voltage_map = voltage_map.sort_values(by='Amplitude (V)', ascending=False)
    # voltage_map.reset_index(drop=True, inplace=True)
    # to always have the voltages in descending order (as used by the following analysis)
    return voltage_map
   
def load_monodisperse_sizes(directory, num_particles):
    sizes = None
    for filename in os.listdir(directory):
        if filename.endswith('.dat'):
            sizes = pd.read_csv(filename, header = None, names = ["Size (um)"])   
            sizes.loc[(sizes["Size (um)"] == 0)] = np.nan   
            size_mean = sizes["Size (um)"].mean()
            sizes.loc[(sizes["Size (um)"] > size_mean*2)] = np.nan # filter out double particles
            size_mean = sizes["Size (um)"].median() # use median instead of mean to get the center of the distribution
            sizes.loc[(sizes["Size (um)"] > size_mean*1.02)] = np.nan # filtering out 
            sizes.loc[(sizes["Size (um)"] < size_mean*0.98)] = np.nan # erroneously determined sizes
            # sizes.loc[(sizes["Size (um)"] > size_mean*1.2)] = np.nan # filtering out 
            # sizes.loc[(sizes["Size (um)"] < size_mean*0.8)] = np.nan # erroneously determined sizes
            # size_mean = sizes["Size (um)"].mean()
            
            size_mean = sizes["Size (um)"].mean()
            sizes = sizes * Phys_properties['diameter_sample'] / size_mean # um; rescale to correct for difraction-induced error
            break    
    if sizes is None:
        sizes = pd.DataFrame(Phys_properties['diameter_sample'], index=np.arange(num_particles), columns=["Size (um)"]) 
    return sizes
    
def load_polydisperse_sizes(directory):
    sizes = None
    for filename in os.listdir(directory):
        if filename.endswith('.dat'):
            sizes = pd.read_csv(filename, header = None, names = ["Size (um)"])
            sizes = sizes / Phys_properties['nm_per_px_def'] * Phys_properties['nm_per_px'] * 2 / 1000  # um
            sizes.loc[(sizes["Size (um)"] == 0)] = np.nan 
            return sizes
        
def load_heatmap():
    heatmap = pd.read_excel(dir_calibr+"\Heatmap.xlsx",index_col=[0])
    x = heatmap.index.to_numpy(dtype=float)
    y = heatmap.columns.to_numpy(dtype=float)  
    f = heatmap.values
    heatmap_min = f.min()
    return x, y, f, heatmap_min
        
def load_data(dir_analysis):   
    print(dir_analysis)
    if calibration_method == "heatmap": # load experimental TV(XY)Z data
        Time, Voltage, XY, Z = load_TVXYZ_data(dir_analysis)
        XY = calculate_xy_positions(XY, dir_analysis)
    else:
        Time, Voltage, Z = load_TVZ_data(dir_analysis)
    num_particles = Z.shape[1]
    annotation = 'Tracking files uploaded: %.0f' %(num_particles)
    print(annotation)
    
    voltage_map = find_voltag_spikes(Voltage) #find voltage application regions
    
    if is_size_known == True: #load sizes of the sample particles
        sizes = load_monodisperse_sizes(dir_analysis, num_particles)
    else:
        sizes = load_polydisperse_sizes(dir_analysis) 
    plot_size_distribution_histogram(sizes)
    r_mean = sizes["Size (um)"].mean() / 2
    r_std = sizes["Size (um)"].std() / 2
    N = sizes["Size (um)"].count()
    annotation = 'Mean size: 𝑑 = ' + str(round(r_mean*2,3)) + ' ± ' + str(round(r_std*2,3)) + ' μm (N = ' + str(N) + ')'
    print(annotation)
    
    if calibration_method == "heatmap": # return values
        return Time, XY, Z, voltage_map, sizes, r_mean, r_std
    else:
        return Time, Z, voltage_map, sizes, r_mean, r_std

def plot_size_distribution_histogram(d):
    d = d.stack().dropna(axis=0).to_numpy(dtype=float)
    
    fig = plt.figure(figsize=(5, 4))  
    ax = fig.add_subplot(1, 1, 1)     
    ax.hist(d, bins = 6)
    ax.set(title = "Size distribution histogram",
       xlabel = "size (μm)",
       ylabel = "count")
    plt.tight_layout()
    plt.show()
    return   
        
#%%
"""
Analysis of the data 
"""
def linear_func(x, a, b):
    y = a * x + b
    return y

def sin_func(x, a, b, l): # amplitude, phase [um], wavelangth
    y = a*np.sin(2*math.pi/l*x-b*(2*math.pi/l))
    return y

def sinking_down_velocity_brenner_fit(z, dz, r, prefactor, density_delta):    
    inner_term = r / (z + r + dz)
    lam = (1 - 9 / 8 * (inner_term) + 1 / 2 * (inner_term)**3 - 57 / 100 * (inner_term)**4 + 1 / 5 * (inner_term)**5 + 7 / 200 * (inner_term)**11 - 1 / 25 * (inner_term)**12)**(-1)
    velocity = prefactor * r**2 * density_delta / lam
    return velocity

def calculate_lambda_brenner (r, h, dz = 0): # h - from the bottom to the bottom of the particle; dz - a systematic shift
    inner_term = r / (h + r + dz)
    lambda_brenner = (1 - 9 / 8 * (inner_term) + 1 / 2 * (inner_term)**3 
           - 57 / 100 * (inner_term)**4 + 1 / 5 * (inner_term)**5 
           + 7 / 200 * (inner_term)**11 - 1 / 25 * (inner_term)**12)**(-1)
    return lambda_brenner

def plot_shooting_up_trajectory(t, z, z_min, r):
    t = t / 1000
    z = z + r - z_min
    fig = plt.figure(figsize=(5, 4))  
    ax = fig.add_subplot(1, 1, 1) 
    ax.plot(t, z, label='Experimental data', color='red')
    ax.axhline(y=0, color='black', linewidth=0.5)
    ax.legend()
    ax.set(title = "Shooting up trajectory",
       xlabel = "time (s)",
       ylabel = "z-position (μm)")
    plt.tight_layout()
    plt.show()
    return

def plot_sin_fit(z, f, params):
    z_fit = np.linspace(0, 20, 51)
    
    fig = plt.figure(figsize=(5, 4))  
    ax = fig.add_subplot(1, 1, 1) 
    ax.plot(z, f, label='Experimental data', color='red')
    ax.plot(z_fit, sin_func(z_fit, params[0], params[1], params[2]), color='dimgray', label='Sin fit')
    ax.axhline(y=0, color='black', linewidth=0.5)
    ax.legend()
    ax.set(title = "Force distribution",
       xlabel = "z-position (μm)",
       ylabel = "acoustic force (pN)")
    plt.tight_layout()
    plt.show()
    return

def plot_sinking_down_linear_fit(z, t, fit):
    fig = plt.figure(figsize=(5, 4))  
    ax = fig.add_subplot(1, 1, 1) 
    
    ax.plot(t, z, label='Experimental data', color='red')
    ax.plot(t ,fit(t), color='dimgray', label='Linear fit')
    ax.set(title = "Sinking down fit",
       xlabel = "time (s)",
       ylabel = "z-position (μm)",
       ylim = [-1, 22])
    ax.legend()
    plt.tight_layout()
    plt.show()    
    return

def plot_sinking_down_full_fit(z, t, v_down_ra, z_down_ra, t_down_ra, r_m, prefactor, density_delta, average_over):
    
    v_down_ra_um = v_down_ra * 10**6
    z_down_ra_um = z_down_ra * 10**6
    r_um = r_m * 10**6
    z = z + r_um
    t_down_fit = np.diff(t_down_ra)
    z_down_fit = z_down_ra_um
    z_down_fit[0] = z_down_ra_um[0] + r_um
    for i in range(1, len(v_down_ra_um)):
        z_down_fit[i] = z_down_fit[i-1] + t_down_fit[i-1]*v_down_ra_um[i]
        
    fig = plt.figure(figsize=(5, 4))  
    ax1 = fig.add_subplot(1, 1, 1) 
    ax2 = ax1.twinx()
    ax1.plot(t, z, color='red',label='Sinking down trajectory')
    ax1.plot(t_down_ra, z_down_fit, color='darkred',label='Sinking down fit')
    ax2.plot(t_down_ra, v_down_ra_um, color='lightgrey', alpha = 1 ,label='Velocity (av. over '+str(average_over)+' points)')
    ax2.plot(t_down_ra, sinking_down_velocity_brenner_fit(z_down_ra, 0, r_m, prefactor, density_delta)*1000000,  color='dimgrey', label = 'Velocity fit')
        
    ax1.set(title = "Sinking down fit",
       xlabel = "time (s)",
       ylabel = "z-position (μm)",
       ylim = [-1, 22])
    ax2.set(ylabel = "𝑣 sinking down (μm/s)")
    ax1.legend()
    plt.tight_layout()
    plt.show()  
    return

def plot_sinking_down_forces_fit(z, t, dz, r_m, prefactor, density_delta, \
        sinking_down_forces_fit_model):
    r_um = r_m * 10**6     
      
    fig = plt.figure(figsize=(5, 4))  
    ax = fig.add_subplot(1, 1, 1) 
    ax.plot(t, z + r_um, color='red',label='Sinking down trajectory')
    ax.plot(t, sinking_down_forces_brenner_fit(t, dz, r_m, prefactor, density_delta, z[0]) + r_um, color='dimgrey',label='Sinking down fit')
    # ax2.plot(t_down_ra, v_down_ra_um, color='lightgrey', alpha = 1 ,label='Velocity (av. over '+str(average_over)+' points)')
    # ax2.plot(t_down_ra, sinking_down_velocity_brenner_fit(z_down_ra, 0, r_m, prefactor, density_delta)*1000000,  color='dimgrey', label = 'Velocity fit')
        
    ax.set(title = "Sinking down fit",
       xlabel = "time (s)",
       ylabel = "z-position (μm)",
       ylim = [-1, 22])
    ax.legend()
    plt.tight_layout()
    plt.show()  
    return

"""
???? If the z_up at the bottom is WRONG (not 0), subtratct absolute min (Z.min()) instead of the local one (z_up_start = Z.iloc[start-5:start+10,i].min()) 
"""
def find_xy_up_values(start, end, i):
    if calibration_method == "heatmap":
        x_up = XY.iloc[start:end, i * 2].to_numpy()
        y_up = XY.iloc[start:end, i * 2 + 1].to_numpy()
    else:
        x_up = 'Nan'
        y_up = 'Nan'
    return x_up, y_up
    
def find_voltage_normalized_force(x_up, y_up, z_up, z_up_start, time_up, h_jump, r_um, voltage):
    
    #setting up z-interval for fitting
    if calibration_method == "heatmap":
        lower_limit = 0.10
        upper_limit = 0.90
        min_num_points = 15
    else:
        lower_limit = 0.3
        upper_limit = 0.7
        min_num_points = 4
        
    z_up_fit_start = h_jump * lower_limit
    z_up_fit_end = h_jump * upper_limit
    z_up = z_up - z_up_start
    
    # finding the starting point index for the fit
    start_fit = 0
    for i in range(len(z_up)):
        if z_up[i] > z_up_fit_start:
            start_fit = i - 1 # the last point below the lower_limit
            break
        
    # finding the ending point index for the fit
    end_fit = 0
    for i in range(len(z_up)):
        if z_up[i] > z_up_fit_end:
            end_fit = i # the first point above the upper_limit
            break
    
    if end_fit - start_fit < min_num_points or start_fit == -1:        
        force_acoustic = np.nan
        x_up_position = np.nan
        y_up_position = np.nan
        return force_acoustic, x_up_position, y_up_position
    z_up_for_fit = z_up[start_fit:end_fit] # um
    time_up_for_fit = time_up[start_fit:end_fit] / 1000 # s
    r_m = r_um / 1000000 # m
    if calibration_method == "heatmap":
        if (    any(x_up < 0) # no xy-tracking errors
            or any(y_up < 0)
            or any(x_up > Phys_properties['FOV_x_um'])
            or any(y_up > Phys_properties['FOV_y_um'])  ):
            
            force_acoustic = np.nan
            x_up_position = np.nan # filter out loss ofxy-tracking
            y_up_position = np.nan
        
        else:
            velocity_up = calculate_velocity(z_up_for_fit, time_up_for_fit) # um/s

            z_up_for_fit = z_up_for_fit[1:-1] # m, at the same positions as velocity_up
            force_acoustic_array = convert_velocity_to_force(velocity_up, r_m, z_up_for_fit) # pN
            
            z_up_for_fit = z_up_for_fit + r_um # to have the bottom as z = 0
            h = h_jump + r_um
            sin_fit_values, residues = shooting_up_sin_fit(z_up_for_fit, force_acoustic_array, h)
    
            # plot_sin_fit(z_up_for_fit, force_acoustic_array, sin_fit_values)

            force_acoustic = sin_func(8, sin_fit_values[0], sin_fit_values[1], sin_fit_values[2]) # at the bottom (z = 0)        
            
            idx_mid = start_fit + (end_fit - start_fit) // 2
            x_up_position = x_up[0] # at the starting position (i.e. at the bottom)
            y_up_position = y_up[0]
    
    else:               
        velocity_up = shooting_up_linear_fit(z_up_for_fit, time_up_for_fit) # um/s  
        z_up_mid = h_jump * 0.5 # um
        force_acoustic = convert_velocity_to_force(velocity_up, r_m, z_up_mid)
        x_up_position = np.nan # x and y are not used for "average" response method
        y_up_position = np.nan
    
    if is_size_known == True:
            force_acoustic = force_acoustic * (Phys_properties['radius_sample'] / r_um)**3 # size correction for homogeneous samples

    force_acoustic_normalized = force_acoustic / voltage**2 # pN/V^2
    return force_acoustic_normalized, x_up_position, y_up_position
        
"""
???? How to calculate the force at the same height for linear fit
"""    

def calculate_velocity(z, t):
    dz = z[2:] - z[:-2]
    dt = t[2:] - t[:-2]
    velocity_up = dz / dt 
    return velocity_up
    
def shooting_up_sin_fit(z_up_for_fit, force_acoustic, h):
    velocity_sound = Phys_properties['velocity_sound'] # m/s
    frequency_resonance = Phys_properties['frequency_resonance'] # MHz
    wavelength = velocity_sound / frequency_resonance / 2 # um
    phase = h - wavelength / 2 # um
    fitfunc = sin_func
    p, pcov = curve_fit(fitfunc, z_up_for_fit, force_acoustic,
                         p0=[1.0, phase, wavelength],
                         bounds=((0     , phase-0.001*wavelength, wavelength*0.995),
                                 (np.inf, phase+0.001*wavelength, wavelength*1.005)))
    return p, pcov

def shooting_up_linear_fit(z_up_for_fit, time_up_for_fit):
    coef, cov = curve_fit(linear_func,
                    time_up_for_fit, z_up_for_fit, 
                    bounds=([0, -np.inf], [np.inf, np.inf]))
    velocity_up = np.poly1d(coef)[1] # single value    
    return velocity_up

def convert_velocity_to_force(velocity_up, r_m, z_up_mid):
    velocity_up = velocity_up / 1000000 # m/s
    z_up_mid = z_up_mid / 1000000 # m
    force_gravity = 4 / 3 * math.pi * (r_m)**3 * Phys_properties['g'] * Phys_properties['density_sample'] # N
    force_buoyancy = - 4 / 3 * math.pi * (r_m)**3 * Phys_properties['g'] * Phys_properties['density_medium'] # N
    lambda_brenner = calculate_lambda_brenner(r_m, z_up_mid) # dimensionless
    # print(lambda_brenner)
    force_stokes_prefactor = 6 * math.pi * r_m * Phys_properties['viscosity_medium'] * lambda_brenner # N / (m/s)
    force_acoustic = (velocity_up * force_stokes_prefactor + force_gravity + force_buoyancy) # N
    force_acoustic = force_acoustic * 10**12 # pN
    return force_acoustic

def find_fitting_interval(z, lower_limit, upper_limit):
    z_max = z.max()
    z_down_fit_start = z_max * upper_limit
    z_down_fit_end = z_max * lower_limit    
    
    
    # finding the starting point index for the fit
    start_fit = 0
    for i in range(1, len(z)):
        if z[i] < z_down_fit_start:
            start_fit = i # the first point below the upper_limit
            break
        
    # finding the ending point index for the fit
    end_fit = 0
    for i in range(len(z)):
        if z[i] < z_down_fit_end:
            end_fit = i # the first point below the lower_limit
            break
    
    return start_fit, end_fit

def find_density_from_v_sinking_down(v, r_um):
    r_m = r_um / 10**6 # m
    v_ms = v / 10**6 # m/s
    density = - 9 * Phys_properties['viscosity_medium'] * v_ms / \
        (2 * Phys_properties['g'] * r_m**2) + Phys_properties['density_medium']
    return density
 
def sinking_down_linear_fit(z, t, start_fit, end_fit, r_um):               
    r_m = r_um / 10**6 # m
    z_mid = (z[start_fit] + z[end_fit]) / 2 # um
    z_mid = z_mid / 10**6 # m
    lambda_brenner = calculate_lambda_brenner(r_m, z_mid)
        
    z = z + r_um
    z_for_fit = z[start_fit:end_fit]
    t_for_fit = t[start_fit:end_fit]
    
    coef, cov = curve_fit(linear_func, t_for_fit, z_for_fit, bounds=([-np.inf,-np.inf], [np.inf, np.inf]))
    v_sinking_down = np.poly1d(coef)[1] * lambda_brenner
    
    # uncomment to plot the fitted sinking down trajectory
    # fit = np.poly1d(coef)     
    # plot_sinking_down_linear_fit(z, t, fit)
        
    return v_sinking_down

def sinking_down_full_fit(z, t, start_fit, end_fit, r_um):
    r_m = r_um / 10**6 # m
    prefactor = - 2 / 9 * Phys_properties['g'] / Phys_properties['viscosity_medium']
    
    z_range = z[start_fit:end_fit]
    t_range = t[start_fit:end_fit]
    
    v_down = np.diff(z_range)/np.diff(t_range) # um/s
    t_down = (t_range[1:]+t_range[:-1])/2 # s
    z_down = (z_range[1:]+z_range[:-1])/2 # um
    
    average_over = math.ceil((end_fit - start_fit) / 10)
    v_down_ra = np.convolve(v_down, np.ones(average_over), 'valid') / average_over #rolling average in um/s
    t_down_ra = np.convolve(t_down, np.ones(average_over), 'valid') / average_over #rolling average in s                        
    z_down_ra = np.convolve(z_down, np.ones(average_over), 'valid') / average_over #rolling average in um
    
    v_down_ra = v_down_ra / 10**6 # m/s
    z_down_ra = z_down_ra / 10**6 # m
    
    sinking_down_full_fit_model = Model(sinking_down_velocity_brenner_fit,independent_vars=['z'])
    
    result = sinking_down_full_fit_model.fit(v_down_ra, z=z_down_ra, dz=Parameter('dz', value=0, vary=False), \
        r=Parameter('r', value=r_m, vary=False), prefactor=Parameter('prefactor', value=prefactor, vary=False), density_delta= 100)
    
    density_delta = result.params['density_delta'].value
    density = density_delta + Phys_properties['density_medium']
    
    # uncomment to plot the fitted sinking down trajectory and velocity fit
    # plot_sinking_down_full_fit(z, t, v_down_ra, z_down_ra, t_down_ra, r_m, \
    #     prefactor, density_delta, average_over)
    
    return density

def sinking_down_forces_brenner_fit(t, dz, r, prefactor, density_delta, z_start):
    z_fit = np.arange(len(t),dtype=float) # m
    v_fit = np.arange(len(t)-1,dtype=float) # m/s
    z_fit_um = np.arange(len(t),dtype=float) # um
    z_fit_um[0] = z_start
    z_fit[0] = z_fit_um[0] / 10**6
    
    for i in range(1, len(z_fit)):
        v_fit[i-1] = prefactor * r**2 * density_delta / calculate_lambda_brenner(r, z_fit[i-1],dz)
        z_guess = z_fit[i-1] + v_fit[i-1] * (t[i] - t[i-1])
        v_fit[i-1] = prefactor * r**2 * density_delta / calculate_lambda_brenner(r, z_guess,dz)
        z_fit[i] = z_fit[i-1] + v_fit[i-1] * (t[i] - t[i-1])
    
    z_fit_um = z_fit * 10**6
    
    return z_fit_um   

def sinking_down_forces_fit(z, t, start_fit, end_fit, r_um):
    r_m = r_um / 10**6 # m
    prefactor = - 2 / 9 * Phys_properties['g'] / Phys_properties['viscosity_medium']
    dz_shift = -0.2 / 10**6 # m
    
    z_range = z[start_fit:end_fit]
    t_range = t[start_fit:end_fit]
    
    sinking_down_forces_fit_model = Model(sinking_down_forces_brenner_fit, independent_vars=['t'])
    
    density = np.nan
    
    try:
        result = sinking_down_forces_fit_model.fit(z_range, t = t_range, dz=Parameter('dz', value=dz_shift, vary=False), \
            r=Parameter('r', value=r_m, vary=False), prefactor=Parameter('prefactor', value=prefactor, vary=False), \
            z_start=Parameter('z_start', value=z_range[0], vary=False), density_delta= 10)
        density_delta = result.params['density_delta'].value
        density = density_delta + Phys_properties['density_medium']
    except:
        pass
    
    
    # uncomment to plot the fitted sinking down trajectory and velocity fit
    # plot_sinking_down_forces_fit(z, t, dz_shift, r_m, prefactor,\
    #     density_delta, sinking_down_forces_fit_model)
    
    return density

def calculate_density(z, h, t, r_um):
    z = z - z.min()
    t = t / 1000 # s
    if density_method == "linear":
        upper_limit = 0.95
        lower_limit = 0.6
        start_fit, end_fit = find_fitting_interval(z, lower_limit, upper_limit)
        if (end_fit - start_fit) < 1:
            return np.nan
        v_sinking_down = sinking_down_linear_fit(z, t, start_fit, end_fit, r_um) # um/s
        density = find_density_from_v_sinking_down(v_sinking_down, r_um)
    
    elif density_method == "full": 
        upper_limit = 0.95
        lower_limit = 0.05
        start_fit, end_fit = find_fitting_interval(z, lower_limit, upper_limit)
        density = sinking_down_full_fit(z, t, start_fit, end_fit, r_um)
        
    elif density_method == "forces":
        upper_limit = 0.95
        lower_limit = 0.05
        start_fit, end_fit = find_fitting_interval(z, lower_limit, upper_limit)
        density = sinking_down_forces_fit(z, t, start_fit, end_fit, r_um)
           
    return density        

def analyze_data():
    quality_up = pd.DataFrame(index=range(num_voltages),columns=range(num_particles)) 
    # is the shooting up trajectory suitable for the analysis
    quality_down = pd.DataFrame(index=range(num_voltages),columns=range(num_particles))  
    #  is the sinking down trajectory suitable for the analysis
    z_up_min = pd.DataFrame(index=range(num_voltages),columns=range(num_particles))  
    # min z during a voltage application
    z_up_max = pd.DataFrame(index=range(num_voltages),columns=range(num_particles))  
    # max z during a voltage application
    force_values = pd.DataFrame(index=range(num_voltages),columns=range(num_particles))  
    # force values stored to later generate a heatmap
    densities = pd.DataFrame(index=range(num_voltages),columns=range(num_particles))  
    # densities calculated from each sinking down event
    x_up_values = pd.DataFrame(index=range(num_voltages),columns=range(num_particles))  
    # x-position values storred to later generate a heatmap
    y_up_values = pd.DataFrame(index=range(num_voltages),columns=range(num_particles))  
    # y-position values storred to later generate a heatmap
    
    for i in range(num_particles):
    # for i in [120]:    
        # print('Particle ' + str(i))
        # for j in range(2,3):
        for j in range(num_voltages):
            voltage = voltage_map.loc[j]['Amplitude (V)']
            start = int(voltage_map.loc[j]['Start'])
            end = int(voltage_map.loc[j]['End'])
            z_up_min.at[j,i] = Z.iloc[0:10,i].mean()
            # z_up_min.at[j,i] = Z.iloc[start-5:start+5,i].mean()
            z_up_max.at[j,i] = Z.iloc[start:end,i].max()
            
            z_up = Z.iloc[start:end,i].to_numpy() # z-trajectory throughout the voltage application
            dz_up = np.diff(z_up, n=1) # difference between the adjucent z points
            time_up = Time.iloc[start:end].to_numpy() # time throughout the voltage application
            time_up = time_up - time_up.min()
            x_up, y_up = find_xy_up_values(start, end, i)
            
            r_um = sizes.loc[i]['Size (um)'] / 2 # um
            h_jump = z_up_max.iloc[j,i] - z_up_min.iloc[j,i] # um
            h_expected = Phys_properties['acoustic_node_position'] - r_um # um
            
            # plot_shooting_up_trajectory(time_up, z_up, z_up_min.at[j,i], r_um)
            
            
            if ( pd.isnull(sizes.iloc[i,0]) # the size of the particle is not knonw
                    or any(z_up == 0) or any(z_up == -40) # extreme values
                    or any(z_up < z_up_min.iloc[j,i] - 1) # unexpected motion downward
                    or h_jump < h_expected * 0.7 # insufficient shooting up height
                    or h_jump > h_expected * 2 # overshooting up
                    or any(abs(dz_up) > 7)  ): # errors of tracking / very fast motion                  
                quality_up.at[j,i] = 0
                quality_down.at[j,i] = 0
            else:
                quality_up.at[j,i] = 1
                force_values.at[j,i], x_up_values.at[j,i], y_up_values.at[j,i] = \
                        find_voltage_normalized_force(x_up, y_up, z_up, 
                        z_up_min.iloc[j,i], time_up, h_jump, r_um, voltage)
                
                start = int(voltage_map.loc[j]['End']) # finding the begining of the sinking down trajectory
                if j != num_voltages - 1: # finding the end of the sinking down trajectory
                    end = int(voltage_map.loc[j+1]['Start'])
                else:
                    end = int(Z.last_valid_index())
                time_down = Time.iloc[start:end].to_numpy() 
                # print(start, end)
                time_down = time_down - time_down.min()  

                z_down = Z.iloc[start:end,i].to_numpy()
                dz_down = np.diff(z_down, n=1) # difference between the adjucent z points
                
                if (    any(z_up == 0) or any(z_up == -40) # extreme values
                        or z_down.max() - z_down.min() < h_expected - r_um # insufficient sinking down depth
                        or any(abs(dz_down) > 5)    ): # errors of tracking / very fast motion
                    quality_down.at[j,i] = 0
                else:
                    quality_down.at[j,i] = 1
                    densities.at[j,i] = calculate_density(z_down, h_jump, time_down, r_um)     

    return force_values, x_up_values, y_up_values, densities

#%%
def recalculate_full_acoustic_force(f, v):
    v_squared = v['Amplitude (V)']**2    
    f_ac = f.multiply(v_squared, axis="index")
    return f_ac, v_squared

def filter_out_particles(f, n):
    threshold = n
    columns = f.count()
    columns = columns[columns >= threshold].index.tolist()
    return columns
    # only leave the particles with >= n shooting up events

def plot_density_distribution_histogram(rho) :
    
    
    rho = rho.dropna().to_numpy(dtype=float)
    mu, sigma = scipy.stats.norm.fit(rho)
    x = np.linspace(1055, 1070, num=21)
    best_fit_line = scipy.stats.norm.pdf(x, mu, sigma)

    print(mu, sigma)
    # rho = y.mean()
    # print(rho)
    fig = plt.figure(figsize=(5, 4))  
    ax = fig.add_subplot(1, 1, 1)     
    ax.hist(rho, bins = 5)
    ax.plot(x, best_fit_line*12*4.5, color = 'red')
    ax.set(title = "Density distribution histogram",
       xlabel = "Density (μm)",
       ylabel = "count",
       xlim = [1055, 1070])
    plt.tight_layout()
    plt.show()

    # fig = plt.figure(figsize=(5, 4))  
    # ax = fig.add_subplot(1, 1, 1)     
    # ax.hist(rho, bins = 6)
    # ax.set(title = "Density distribution histogram",
    #    xlabel = "Density (kg/m³)",
    #    ylabel = "count")
    # plt.tight_layout()
    # plt.show()
    return

def plot_individual_force_vs_voltage_squared(force_values, voltage_map, y_label):
    columns = filter_out_particles(force_values, 2)
    force_filtered = force_values.loc[:,columns]
    f_ac, v_squared = recalculate_full_acoustic_force(force_filtered, voltage_map) # pN not voltage-squared normalized 
    num_lines = len(columns)
    cm_subsection = np.linspace(0.0, 1.0, num_lines) 
    colors = [ cm.viridis(x) for x in cm_subsection ]
    fig = plt.figure(figsize=(5, 4))  
    ax = fig.add_subplot(1, 1, 1) 
    for i, f_idx in enumerate(f_ac): 
        # ax.plot(v_squared, f_ac[f_idx], color = colors[i], marker='o', linestyle='none')
        # coef, cov = fit_average_force_vs_voltage_squared(v_squared, f_ac[f_idx], 8, 6)
        # ax.plot(v_squared, linear_func(v_squared, coef[0], coef[1]), color='red', label='Linear fit')
        
        coef = 0
        count = 0
        # if f_ac[f_idx][8] != 0:
        #     count = count + 1
        #     coef = coef + f_ac[f_idx][8] / v_squared[8]
        if f_ac[f_idx][7] != 0:
            count = count + 1
            coef = coef + f_ac[f_idx][7] / v_squared[7]
        if f_ac[f_idx][6] != 0:
            count = count + 1
            coef = coef + f_ac[f_idx][6] / v_squared[6]
        if f_ac[f_idx][5] != 0:
            count = count + 1
            coef = coef + f_ac[f_idx][5] / v_squared[5]
        ax.plot(v_squared, linear_func(v_squared, coef/count, 0), color='grey', label='Linear fit')
        ax.plot(v_squared, f_ac[f_idx], marker='o', linestyle='none', color = colors[i])

    ax.set(title = "Individual responses",
       xlabel = "voltage squared (V²)",
       ylabel = y_label)
       # xlim = [-1, 12.5],
       # ylim = [-0.5, 5])
    plt.tight_layout()
    plt.show()
    return

def plot_average_force_vs_voltage_squared(v_squared, f, f_std, coef, y_label):
    # print(v_squared, f, coef)
    fig = plt.figure(figsize=(5, 4))  
    ax = fig.add_subplot(1, 1, 1) 
    ax.errorbar(v_squared, f, f_std, label='Experimental data', color='red', marker='o')
    ax.plot(v_squared, linear_func(v_squared, coef[0], coef[1]), color='dimgray', label='Linear fit')
    ax.legend()
    ax.set(title = "Average response",
       xlabel = "voltage squared (V²)",
       ylabel = y_label,
       xlim = [-1, 12.5],
       ylim = [-1, 10])
    plt.tight_layout()
    plt.show()
    return

def plot_average_and_individual_force_vs_voltage_squared(f, coef, v_squared, y_label):
    columns = filter_out_particles(f, 2)
    force_filtered = f.iloc[:,columns]
    num_lines = len(columns)

    fig = plt.figure(figsize=(5, 4))  
    ax = fig.add_subplot(1, 1, 1) 
    for i, f_idx in enumerate(force_filtered): 
        ax.plot(v_squared, force_filtered[f_idx], color = 'grey')
    ax.plot(v_squared, linear_func(v_squared, coef[0], coef[1]), color='red', label='Linear fit')
    ax.set(title = "Individual responses",
       xlabel = "voltage squared (V²)",
       ylabel = y_label,
       xlim = [-1, 12.5],
       ylim = [-0.5, 5])
    plt.tight_layout()
    plt.show()
    return

def plot_error_of_the_fit(f_ac_ratio_average, v_squared, f_ac_ratio_average_fit):
    errors = (f_ac_ratio_average - v_squared * f_ac_ratio_average_fit) / f_ac_ratio_average
    fig = plt.figure(figsize=(5, 4))  
    ax = fig.add_subplot(1, 1, 1) 
    ax.plot(v_squared, errors, color='red', marker='o', linestyle='none')
    plt.tight_layout()
    plt.show()
    return

def plot_binned_individual_response(f_ac, f_ac_std, size_3):
    x = size_3
    y = f_ac
    print(y
          )
    bins_list = [0, 63, 102, 130 , 260, 400]
    
    n, _ = np.histogram(x, bins=bins_list)
    sy, _ = np.histogram(x, bins=bins_list, weights=y)
    sy2, _ = np.histogram(x, bins=bins_list, weights=y*y)
    sx, _ = np.histogram(x, bins=bins_list, weights=x)
    sx2, _ = np.histogram(x, bins=bins_list, weights=x*x)
    mean_y = sy / n
    mean_x = sx / n
    std_y = np.sqrt(sy2/n - mean_y*mean_y)
    std_x = np.sqrt(sx2/n - mean_x*mean_x)
    
    coef, cov = curve_fit(linear_func,mean_x,mean_y, bounds=([-np.inf, -0.01], [np.inf, 0.01]))
    slope = coef[0]
    slope_err = np.sqrt(np.diag(cov))[0]
    fig = plt.figure(figsize=(5, 4))  
    ax = fig.add_subplot(1, 1, 1) 
    ax.plot(x, y, 'bo', color='lightgrey',zorder=0)
    ax.errorbar(mean_x, mean_y, xerr = std_x, yerr=std_y, fmt='o', linestyle='none')
    ax.plot([0,x.max()*1.1],[0,x.max()*1.1*slope],color='steelblue')
    ax.fill_between([0,x.max()*1.1], [0,x.max()*1.1*(slope-np.sqrt(np.diag(cov))[0])], [0,x.max()*1.1*(slope+np.sqrt(np.diag(cov))[0])], color='steelblue', alpha=0.15)
    
    ax.set(title = "Individual binned responses",
       xlabel = "size cubed (μm³)",
       ylabel = "force ratio")
    plt.tight_layout()
    plt.show()
    
    size_norm_factor = (Phys_properties['radius_calibr']*2)**3
    slope_size_norm = slope * size_norm_factor
    slope_size_norm_err = slope_err * size_norm_factor
    contrast_factor_sample, contrast_factor_sample_std = \
            find_contrast_factor(slope_size_norm, slope_size_norm_err, Phys_properties['radius_calibr'], 0) # *4.65**3 to make the slope dimensionless (4.65, because it was the calibration)
    compr_sample, compr_sample_std = find_compr(contrast_factor_sample, contrast_factor_sample_std)
    annotation = 'Contrast factor: 𝜑 = ' + str(round(contrast_factor_sample,3)) + ' ± ' + str(round(contrast_factor_sample_std,3))
    print(annotation)
    
    annotation = 'Compressibility: 𝛽 = ' + str(round(compr_sample,3)) + ' ± ' + str(round(compr_sample_std,3)) + ' 10⁻¹⁰ Pa⁻¹'
    print(annotation)
    return
    

def plot_norm_force_histogram(f):

    fig = plt.figure(figsize=(5, 4))  
    ax = fig.add_subplot(1, 1, 1)     
    ax.hist(f, bins = 10)
    ax.set(title = "Force distribution histogram",
       xlabel = "normalized force (pN/V²)",
       ylabel = "count")
    plt.tight_layout()
    plt.show()
    return

def plot_compressibility_distribution(compr_sample):
    fig = plt.figure(figsize=(5, 4))  
    ax = fig.add_subplot(1, 1, 1)
    ort="h"
    sigma = .4
    ax = pt.RainCloud(x = compr_sample, bw = sigma, width_viol = .3, ax = ax, orient = ort, move = 0.25, alpha = 0.5)
    ax.set(title = "Compressibility distribution",
       ylabel = "compressibility (10⁻¹⁰ Pa⁻¹)",
       ylim = [0,5])
    plt.tight_layout()
    plt.show()
    return

def plot_TRAP(x_TRAP, x_TRAP_av,x_TRAP_av_std, v_times, title, y_label): # x - force, compressibility, etc.
    fig = plt.figure(figsize=(5, 4))  
    ax = fig.add_subplot(1, 1, 1) 
    for i, x_TRAP_idx in enumerate(x_TRAP):
        ax.plot(v_times, x_TRAP[x_TRAP_idx], color = 'grey')
    ax.plot(v_times, x_TRAP_av, color='red', label='Linear fit')
    ax.fill_between(v_times, x_TRAP_av-x_TRAP_av_std, x_TRAP_av+x_TRAP_av_std, alpha=0.05, color = 'red')
    print(x_TRAP_av_std)
    ax.set(title = title,
       xlabel = "time (min)",
       ylabel = y_label)
    plt.tight_layout()
    plt.show()
    return

def smoothen_heatmap(f, x, y):
    for i in range(len(f)-1):
        nearby = [f[i]]
        nearby_idx = [i]
        for j in range(i+1,len(f)):
            distance = np.sqrt((x[i] - x[j])**2 + (y[i] - y[j])**2)
            if distance < 30:
                nearby = np.append(nearby,f[j])
                nearby_idx = np.append(nearby_idx,j)
        if len(nearby_idx) > 1:
            nearby_mean_value = np.mean(nearby)
            for idx in nearby_idx:
                f[idx] = nearby_mean_value
    return f
    
def produce_heatmap(f):
  
    f_heatmap = f.iloc[:,0].to_numpy(dtype=float)
    x_heatmap = f.iloc[:,1].to_numpy(dtype=float)
    y_heatmap = f.iloc[:,2].to_numpy(dtype=float)

    # f_heatmap_av = np.mean(f_heatmap)
    # f_heatmap = (f_heatmap - f_heatmap_av)/f_heatmap_av
    # print(f_heatmap_av)
    # print(f_heatmap)
    
    f_heatmap_smooth = smoothen_heatmap(f_heatmap, x_heatmap, y_heatmap)
    
    FOV_x = Phys_properties['FOV_x_um']
    FOV_y = Phys_properties['FOV_y_um']
    f_min = f_heatmap_smooth.min()
    f_max = f_heatmap_smooth.max()
    x_grid, y_grid = np.mgrid[0:FOV_x:126j, 0:FOV_y:101j]
    heatmap = griddata((x_heatmap, y_heatmap), f_heatmap_smooth, (x_grid,y_grid), method='linear',fill_value=0)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    plt.imshow(heatmap.T, extent=[0,FOV_x,0,FOV_y], origin='lower',cmap='rainbow')
    # ax.scatter(x_heatmap, y_heatmap, s=f_heatmap_smooth*10,facecolors='none',edgecolors='dimgrey')

    ax.set(title = "Field intensity heatmap (pN/V²)", 
            xlabel = "x (μm)", 
            ylabel = "y (μm)",
            xlim = [0,FOV_x],
            ylim = [0,FOV_y])
    
    # # uncomment to add particle indices    
    # for i, particle_num in enumerate(combined_fxy_heatmap.index):
    #     ax.text(x_heatmap[i] + 5, y_heatmap[i] + 1, str(combined_fxy_heatmap.iloc[i,3]) + '-' + str(particle_num), fontsize=7)
        
    plt.clim(0,1.5)
    plt.colorbar(fraction=0.03777, pad=0.02)
    plt.show()
    
    plot_norm_force_histogram(f_heatmap_smooth) # plot distribution of forces

    annotation = 'Heatmap statistics: N = ' + str(f_heatmap_smooth.size)
    print(annotation)
    col_idx = pd.Series(y_grid[1,:]) # export the heatmap
    row_idx = pd.Series(x_grid[:,1]) # x - rows, y - columns
    heatmap_out = pd.DataFrame(heatmap, index=row_idx, columns=col_idx, dtype=None, copy=False)
    heatmap_out.to_excel('Heatmap.xlsx')         
    print('Heatmap mean f = ', heatmap_out.stack().mean(), '; std f = ', heatmap_out.stack().std())
    r_calibr = r_mean
    r_calibr_std = r_std
    contr_factor_calibr = calculate_contr_factor_PS()
    calibration_values = pd.DataFrame([[r_calibr, r_calibr_std, contr_factor_calibr]],)
    file_calibration = "Calibration values.xlsx"
    calibration_values.to_excel(file_calibration, index=False)
    return

def find_optimal_fit_cut_offs(f, v_squared, n, l ):
    # find the optimal cut offs by minimizing the error of the intersection with the y-axis
    # n+1 - min num of points; l - range of the lower cut-off variation
    # coef[0] - slope, coef[1] - intersection
    threshold = n # min number of points
    lower_range = l # these many lowest voltages are tried as cut offs
    higher_range = num_voltages - lower_range - threshold + 1 # same for the highest voltages
    fit_err = [[np.inf] * higher_range for i in range(lower_range)] # rows for lower, cols for highwe
    for i in range(lower_range):
        for j in range(higher_range):
            idx_low = num_voltages - i
            idx_high = j
            try:            
                coef, cov = curve_fit(linear_func,
                    v_squared[idx_high:idx_low], f[idx_high:idx_low], 
                    bounds=([-np.inf, -0.1], [np.inf, 0.1]))
                fit_err[i][j] = np.sqrt(np.diag(cov))[1]
            except:
                pass  
    are_non_inf = (np.min(fit_err) < np.inf)
    if are_non_inf == True:
        idx_min = np.argwhere(fit_err == np.min(fit_err))
        idx_low = num_voltages - idx_min[0,0]
        idx_high = idx_min[0,1] 
        # print(fit_err)
        return idx_low, idx_high
    else:
        return float("NaN"), float("NaN")

def find_average_response(f, v):
    columns = filter_out_particles(f, 2)
    N = len(columns)
    f_filtered = f.iloc[:,columns]
    f_average = f_filtered.mean(axis=1)
    f_average_std = f_filtered.std(axis=1)
    return f_average, f_average_std, N

def fit_average_force_vs_voltage_squared(v_squared, f_ac_average, idx_low, idx_high):
    coef, cov = curve_fit(linear_func,
                    v_squared[idx_high:idx_low], f_ac_average[idx_high:idx_low], 
                    bounds=([-np.inf, -0.1], [np.inf, 0.1]))
    return coef, cov
    
def calculate_contr_factor_PS():
    density_term = (5 * Phys_properties['density_sample'] - 2 * Phys_properties['density_medium']) / \
        (2 * Phys_properties['density_sample'] + Phys_properties['density_medium'])       
    compr_term =  Phys_properties['compr_PS'] / Phys_properties['compr_medium'] 
    contr_factor = density_term - compr_term
    return contr_factor

def average_response_calibration(f, v):
    f_average, f_average_std, N = find_average_response(f, v)
    f_ac_average, v_squared = recalculate_full_acoustic_force(f_average, v) # pN
    f_ac_average_std, v_squared = recalculate_full_acoustic_force(f_average_std, v) # pN
    idx_low, idx_high = find_optimal_fit_cut_offs(f_ac_average, v_squared, 4, 2)
    
    coef, cov = fit_average_force_vs_voltage_squared(v_squared, f_ac_average, idx_low, idx_high)
    
    plot_average_force_vs_voltage_squared(v_squared, f_ac_average, f_ac_average_std, coef, "acoustic force (pN)")

    # uncomment to plot average and individual responses on the same figure
    f_ac, v_squared = recalculate_full_acoustic_force(f, v) # pN
    plot_average_and_individual_force_vs_voltage_squared(f_ac, coef, v_squared, "acoustic force (pN)")
    
    slope_average = coef[0]
    slope_average_std = np.sqrt(np.diag(cov))[0]
    r_calibr = r_mean
    r_calibr_std = r_std
    contr_factor_calibr = calculate_contr_factor_PS()
    calibration_values = pd.DataFrame([[r_calibr, r_calibr_std, contr_factor_calibr,slope_average, slope_average_std]],)
    file_calibration = "Calibration values.xlsx"
    calibration_values.to_excel(file_calibration, index=False)
     
    annotation = 'Force fit ratio: f = ' + str(round(slope_average/slope_average,3)) + ' ± ' + str(round(slope_average_std/slope_average,3)) + ' (N = ' + str(N) + ')'
    print(annotation)

    return

def apply_heatmap_calibration(f, x, y):
    f_ratio = pd.DataFrame(index=range(num_voltages),columns=range(num_particles)) 
    for i in range(num_particles):
        for j in range(num_voltages):
            if pd.notnull(f.iloc[j,i]):
                heatmap_corr = float(heatmap_func([x.iloc[j,i],y.iloc[j,i]]))
                f_ratio.at[j,i] = f.iloc[j,i] / heatmap_corr
                if heatmap_corr < heatmap_min * 1.05:
                    f_ratio.at[j,i] = np.nan
    # how much stronger/weaker the sample particle reacted in comparison to the PS particle of known size

    return f_ratio

def apply_average_calibration(f, f_std):
    f_ratio = f / Phys_properties['calibr']
    f_ratio_std = np.sqrt((f_std / f)**2 + (Phys_properties['calibr_std'] / Phys_properties['calibr'])**2) * f_ratio
    return f_ratio, f_ratio_std

def prepare_analyzed_data_for_TRAP(f, v):
    threshold = num_voltages * 0.6
    columns_total = filter_out_particles(f, threshold) # min number of shooting up events - filter
    columns_first = filter_out_particles(f.iloc[0:1,:], 1) # response at the first voltage - filter
    columns = np.intersect1d(columns_total, columns_first)
    forces_tr = f.iloc[:,columns]
    voltage_times = Time.iloc[(v.loc[:]['Start'] + v.loc[:]['End']) / 2] # ms voltage application mid-time
    voltage_times = voltage_times / 1000 / 60 # min
    return forces_tr, voltage_times, columns

def average_values_TRAP(f):
    N = len(f.columns)
    f_average = f.mean(axis=1)
    f_average_std = f.std(axis=1)
    return f_average, f_average_std, N

def analysis_TRAP(f_ac, voltage_times, r_TRAP):
    f_ac_average, f_ac_average_std, N = average_values_TRAP(f_ac) # average force ratio values at each voltage application
    plot_TRAP(f_ac, f_ac_average, f_ac_average_std, voltage_times, "TRAP: force ratio", "force ratio")
    
    # print(sizes_TRAP)
    f_ac_0 = f_ac / f_ac.iloc[0] # force ratio relative to the first shooting up
    f_ac_average_0, f_ac_average_0_std, N = average_values_TRAP(f_ac_0)
    plot_TRAP(f_ac_0, f_ac_average_0, f_ac_average_0_std, voltage_times, "TRAP: self-normalized force ratio", "self-normalized force ratio")
    
    contrast_factor = (Phys_properties['radius_calibr'] / r_TRAP)**3 * \
        f_ac * Phys_properties['contrast_factor_calibr'] 
    contrast_factor_average, contrast_factor_average_std, N = average_values_TRAP(contrast_factor) 
    plot_TRAP(contrast_factor, contrast_factor_average, contrast_factor_average_std, \
              voltage_times, "TRAP: contrast factor", "contrast factor")
        
    densities_term = (5 * Phys_properties['density_sample'] - 2 * Phys_properties['density_medium']) / \
        (2 * Phys_properties['density_sample'] + Phys_properties['density_medium'])
    compr = (densities_term - contrast_factor) * Phys_properties['compr_medium']
    compr_average, compr_average_std, N = average_values_TRAP(compr) 
    plot_TRAP(compr, compr_average, compr_average_std, \
              voltage_times, "TRAP: compressibility", "compressibility (10⁻¹⁰ Pa⁻¹)")
    print('N = ' + str(N))
    return

def calculate_force_ratio_average_calibr(f, f_std):
    force_ratio_values,force_ratio_values_std = \
        apply_average_calibration(f, f_std) 
        # dimensionless (ratio of pN/V^2 of the sample and the calibration)
    return force_ratio_values, force_ratio_values_std

def find_contrast_factor(force_ratio_values, force_ratio_values_std, r, r_std):
    contrast_factor_sample = (Phys_properties['radius_calibr'] / r)**3 * \
        force_ratio_values * Phys_properties['contrast_factor_calibr'] 
    # contrast_factor_sample_std = np.sqrt((force_ratio_values_std/force_ratio_values)**2 + \
    #     3 * (Phys_properties['radius_calibr_std'] / Phys_properties['radius_calibr'])**2 + \
    #     3 * (r_std / r)**2) * contrast_factor_sample
    
    contrast_factor_sample_std = np.sqrt((force_ratio_values_std/force_ratio_values)**2 + \
        3 * (r_std / r)**2) * contrast_factor_sample
    # contrast_factor_sample_std = np.sqrt((force_ratio_values_std/force_ratio_values)**2) * contrast_factor_sample
    return contrast_factor_sample, contrast_factor_sample_std

def find_compr(contrast_factor_sample, contrast_factor_sample_std):
    densities_term = (5 * Phys_properties['density_sample'] - 2 * Phys_properties['density_medium']) / \
        (2 * Phys_properties['density_sample'] + Phys_properties['density_medium'])
    compr_sample = (densities_term - contrast_factor_sample) * Phys_properties['compr_medium']
    compr_sample_std = np.sqrt((contrast_factor_sample_std/contrast_factor_sample)**2) * compr_sample
    return compr_sample, compr_sample_std

#%%

Phys_properties = set_config_values() #initialize the dictionary with physical properties
folders = select_folders()

if calibration_method == "heatmap":
    heatmap_headers = ['force', 'x (um)', 'y (um)']
    combined_fxy_heatmap = pd.DataFrame(columns=heatmap_headers)

for dir_num, dir_analysis in enumerate(folders):
    os.chdir(dir_analysis)
    print(dir_num, dir_analysis)

    if calibration_method == "heatmap": # depending on the calibration a different set of data is loaded
        Time, XY, Z, voltage_map, sizes, r_mean, r_std = load_data(dir_analysis)
    else:
        Time, Z, voltage_map, sizes, r_mean, r_std = load_data(dir_analysis)
    
    # plot_size_distribution_histogram(sizes)
    num_voltages = len(voltage_map)
    num_particles = Z.shape[1]
    
    if calibration_method == "heatmap" and is_calibration == False: # load previously created heatmap and produce its function
        heatmap_x, heatmap_y, heatmap_f, heatmap_min = load_heatmap() # load the xy matrix of normalized forces
    
        heatmap_func = RegularGridInterpolator((heatmap_x, heatmap_y), heatmap_f, method='linear', fill_value=heatmap_min) 
        #produce a heatmap function
    
    force_values, x_up_values, y_up_values, densities = analyze_data()
    
    columns = filter_out_particles(densities, 1)
    N = len(columns)
    density_mean = densities.mean()
    plot_density_distribution_histogram(density_mean)  
    density_mean.dropna().to_excel("Densities.xlsx")

    density_mean = densities.mean().mean() 
    density_sem = densities.mean().sem() 
    annotation = 'Density: 𝜌 = ' + str(round(density_mean,1)) + ' ± ' + str(round(density_sem,1)) + ' kg/m\u00B3 (N = ' + str(N) + ')' 
    print(annotation)
    if is_density_known == False:
        Phys_properties['density_sample'] = density_mean
    
    # produce calibration
    if is_calibration == True:
        if calibration_method == "heatmap":
            
            f_v_sum = pd.DataFrame()
            f_v = force_values
            f_v.replace(to_replace=np.nan, value=0, inplace = True)
            f_v_sum = f_v.iloc[0,:]+f_v.iloc[1,:]+f_v.iloc[2,:]
            print(f_v_sum[f_v_sum > 0 ].count())

            force_vals = force_values.to_numpy().flatten()
            x_up_vals = x_up_values.to_numpy().flatten()
            y_up_vals = y_up_values.to_numpy().flatten()
            force_vals = pd.DataFrame(force_vals).transpose()
            x_up_vals = pd.DataFrame(x_up_vals).transpose()
            y_up_vals = pd.DataFrame(y_up_vals).transpose()
            fxy_heatmap = [force_vals, x_up_vals, y_up_vals]
            fxy_heatmap = pd.concat(fxy_heatmap, axis=0)
            fxy_heatmap = fxy_heatmap.reset_index(drop=True)
            fxy_heatmap = fxy_heatmap.transpose()
            fxy_heatmap.columns = heatmap_headers
            fxy_heatmap.insert(0, 'Measurement', str(dir_num+1))
            fxy_heatmap.dropna(axis='rows', inplace = True)
            combined_fxy_heatmap = pd.concat([combined_fxy_heatmap,fxy_heatmap])
        if calibration_method == "average": 
            plot_individual_force_vs_voltage_squared(force_values, voltage_map, "acoustic force (pN)")
            average_response_calibration(force_values, voltage_map)
    
    # apply calibration and analyze data
    elif is_calibration == False:
        force_values = force_values * corr_bead
        if is_time_resolved == False:
            if is_size_known == True and calibration_method == "average":     
                # monodisperse average
                plot_individual_force_vs_voltage_squared(force_values, voltage_map, "acoustic force (pN)")
                f_average, f_average_std, N = find_average_response(force_values, voltage_map)
                f_ac_average, v_squared = recalculate_full_acoustic_force(f_average, voltage_map) # pN
                f_average_std, v_squared = recalculate_full_acoustic_force(f_average_std, voltage_map) # pN
                print(f_ac_average, f_average_std)
                idx_low, idx_high = find_optimal_fit_cut_offs(f_ac_average, v_squared, 4, 2)
                coef, cov = fit_average_force_vs_voltage_squared(v_squared, f_ac_average, idx_low, idx_high)
                plot_average_force_vs_voltage_squared(v_squared, f_ac_average, f_average_std, coef, "acoustic force (pN)")
                force_average_fit = coef[0] # pN/V^2
                force_average_fit_std = np.sqrt(np.diag(cov))[0]
                
                force_ratio_values,force_ratio_values_std = \
                    apply_average_calibration(force_average_fit, force_average_fit_std)
                contrast_factor_sample, contrast_factor_sample_std = \
                    find_contrast_factor(force_ratio_values, force_ratio_values_std, r_mean, r_std)
                compr_sample, compr_sample_std = find_compr(contrast_factor_sample, contrast_factor_sample_std)
                
                annotation = 'Force fit ratio: f = ' + str(round(force_ratio_values,3)) + ' ± ' + str(round(force_ratio_values_std,3))
                print(annotation)
                
                # uncomment for Silica/PS measurement, when tthe calibration size was not 4.65 um
                size_corr_force_ratio = force_ratio_values * (Phys_properties['radius_calibr']/2.325)**3
                size_corr_force_ratio_values_std = force_ratio_values_std * (Phys_properties['radius_calibr']/2.325)**3
                annotation = 'Force fit ratio relative to 4.65 um PS beads: f = ' + str(round(size_corr_force_ratio,3)) + ' ± ' + str(round(size_corr_force_ratio_values_std,3))
                print(annotation)
                
                annotation = 'Contrast factor: 𝜑 = ' + str(round(contrast_factor_sample,3)) + ' ± ' + str(round(contrast_factor_sample_std,3)) + ' (N = ' + str(N) + ')'
                print(annotation)
                
                annotation = 'Compressibility: 𝛽 = ' + str(round(compr_sample,3)) + ' ± ' + str(round(compr_sample_std,3)) + ' 10⁻¹⁰ Pa⁻¹ (N = ' + str(N) + ')'
                print(annotation)
            
            elif is_size_known == False and calibration_method == "average":
                # polydisperse average
                f_ac_values, v_squared = recalculate_full_acoustic_force(force_values, voltage_map) # pN
                columns = filter_out_particles(f_ac_values, 1)
                f_ac_filtered = f_ac_values.iloc[:,columns]
                f_ac_fit = pd.DataFrame(index = columns, columns = [['Fit coef (pN/V^2)','Error (pN/V^2)', 'Radius (um)']])
                for i, f_idx in enumerate(f_ac_filtered): 
                    f_ac_for_fit = f_ac_filtered[f_idx]
                    idx_low, idx_high = find_optimal_fit_cut_offs(f_ac_for_fit, v_squared, 0, 3)
                    if np.isfinite(idx_low):
                        coef, cov = fit_average_force_vs_voltage_squared(v_squared, f_ac_for_fit, idx_low, idx_high)
                        f_ac_fit.at[f_idx, 'Fit coef (pN/V^2)'] = coef[0]
                        f_ac_fit.at[f_idx, 'Error (pN/V^2)'] = np.sqrt(np.diag(cov))[0]
                        f_ac_fit.at[f_idx, 'Radius (um)'] = sizes.at[f_idx, 'Size (um)'] / 2
                        
                        # plot_average_force_vs_voltage_squared(v_squared, f_ac_fit, , coef, "acoustic force (pN)")        
                plot_individual_force_vs_voltage_squared(force_values, voltage_map, "acoustic force (pN)")
                
                f_ac_fit.dropna(axis = 0, inplace = True)
                f_ac = f_ac_fit.iloc[:,0].to_numpy(dtype=float)        
                f_ac_std = f_ac_fit.iloc[:,1].to_numpy(dtype=float)
                r_um = f_ac_fit.iloc[:,2].to_numpy(dtype=float)
                
                force_ratio_values,force_ratio_values_std = \
                    apply_average_calibration(f_ac, f_ac_std)        
                contrast_factor_sample, contrast_factor_sample_std = \
                    find_contrast_factor(force_ratio_values, force_ratio_values_std, r_um, 0)
                compr_sample, compr_sample_std = find_compr(contrast_factor_sample, contrast_factor_sample_std)
                
                compr_mean = compr_sample.mean()
                compr_std = compr_sample.std() # /np.sqrt(len(compr_sample)) # sem
                
                N = len(compr_sample)
                annotation = 'Compressibility: 𝛽 = ' + str(round(compr_mean,3)) + ' ± ' + str(round(compr_std,3)) + ' 10⁻¹⁰ Pa⁻¹ (N = ' + str(N) + ')'
                print(annotation)
                
                pd_compr_sample = pd.DataFrame(compr_sample)
                pd_compr_sample.to_excel("Compressibilities.xlsx")
                
                # plot_compressibility_distribution(compr_sample)
                
                # uncomment to plot binned data for PMMA
                # size_3 = (2*r_um)**3
                # print('N = ' + str(len(size_3)))
                # plot_binned_individual_response(force_ratio_values, force_ratio_values_std, size_3)
                
            elif is_size_known == True and calibration_method == "heatmap":     
                # monodisperse heatmap
                force_ratio_values = \
                    apply_heatmap_calibration(force_values, x_up_values, y_up_values)
                    
                plot_individual_force_vs_voltage_squared(force_ratio_values, voltage_map, "Calibrated response (V²)")
                f_ac_ratio, v_squared = recalculate_full_acoustic_force(force_ratio_values, voltage_map) # pN
                f_ac_ratio_average, f_ac_ratio_average_std, N = find_average_response(f_ac_ratio, voltage_map)
                idx_low, idx_high = find_optimal_fit_cut_offs(f_ac_ratio_average, v_squared, 4, 2)
                coef, cov = fit_average_force_vs_voltage_squared(v_squared, f_ac_ratio_average, idx_low, idx_high)
                plot_average_force_vs_voltage_squared(v_squared, f_ac_ratio_average, f_ac_ratio_average_std, coef, "Calibrated response (V²)")
                f_ac_ratio_average_fit = coef[0] # pN/V^2
                f_ac_ratio_average_fit_std = np.sqrt(np.diag(cov))[0]
                
                plot_error_of_the_fit(f_ac_ratio_average, v_squared, f_ac_ratio_average_fit)
                
                annotation = 'Force fit ratio: f = ' + str(round(f_ac_ratio_average_fit,3)) + ' ± ' + str(round(f_ac_ratio_average_fit_std,3))
                print(annotation)
                
                contrast_factor_sample, contrast_factor_sample_std = \
                    find_contrast_factor(f_ac_ratio_average_fit, f_ac_ratio_average_fit_std, r_mean, r_std)
                compr_sample, compr_sample_std = find_compr(contrast_factor_sample, contrast_factor_sample_std)
                
                annotation = 'Compressibility: 𝛽 = ' + str(round(compr_sample,3)) + ' ± ' + str(round(compr_sample_std,3)) + ' 10⁻¹⁰ Pa⁻¹' + ' (N = ' + str(N) + ')'
                print(annotation)
        
            elif is_size_known == False and calibration_method == "heatmap":  
                # polydisperse heatmap  
                force_ratio_values = \
                    apply_heatmap_calibration(force_values, x_up_values, y_up_values)
                    
                f_ac_ratio_values, v_squared = recalculate_full_acoustic_force(force_ratio_values, voltage_map) # pN
                columns = filter_out_particles(f_ac_ratio_values, 3)
                f_ac_ratio_filtered = f_ac_ratio_values.iloc[:,columns]
                f_ac_ratio_fit = pd.DataFrame(index = columns, columns = [['Fit coef (pN/V^2)','Error (pN/V^2)', 'Radius (um)']])
                for i, f_idx in enumerate(f_ac_ratio_filtered): 
                    f_ac_ratio_for_fit = f_ac_ratio_filtered[f_idx]
                    idx_low, idx_high = find_optimal_fit_cut_offs(f_ac_ratio_for_fit, v_squared, 2, 3)
                    if np.isfinite(idx_low):
                        coef, cov = fit_average_force_vs_voltage_squared(v_squared, f_ac_ratio_for_fit, idx_low, idx_high)
                        f_ac_ratio_fit.at[f_idx, 'Fit coef (pN/V^2)'] = coef[0]
                        f_ac_ratio_fit.at[f_idx, 'Error (pN/V^2)'] = np.sqrt(np.diag(cov))[0]
                        f_ac_ratio_fit.at[f_idx, 'Radius (um)'] = sizes.at[f_idx, 'Size (um)'] / 2
                        if coef[0] > 15:
                            f_ac_ratio_fit.at[f_idx, 'Fit coef (pN/V^2)'] = np.nan # discard anomalously high force ratios 
                            # plot_average_force_vs_voltage_squared(v_squared, f_ac_ratio_for_fit, coef, "Calibrated response (V^2)")        
                            
                f_ac_ratio_fit.dropna(axis = 0, inplace = True)
                f_ac_ratio = f_ac_ratio_fit.iloc[:,0].to_numpy(dtype=float)
                f_ac_ratio_std = f_ac_ratio_fit.iloc[:,1].to_numpy(dtype=float)
                r_um = f_ac_ratio_fit.iloc[:,2].to_numpy(dtype=float)
                contrast_factor_sample, contrast_factor_sample_std = \
                    find_contrast_factor(f_ac_ratio, f_ac_ratio_std, r_um, 0)
                compr_sample, compr_sample_std = find_compr(contrast_factor_sample, contrast_factor_sample_std)
        
                compr_mean = compr_sample.mean()
                compr_std = compr_sample.std()
                annotation = 'Compressibility: 𝛽 = ' + str(round(compr_mean,3)) + ' ± ' + str(round(compr_std,3)) + ' 10⁻¹⁰ Pa⁻¹'
                print(annotation)
                
                plot_compressibility_distribution(compr_sample)
        elif is_time_resolved == True:
            if is_size_known == True and calibration_method == "average": 
                f_TRAP, voltage_times, columns = prepare_analyzed_data_for_TRAP(force_values, voltage_map) # filter particles, produce time points array
                f_ac = f_TRAP / Phys_properties['calibr'] # calculate f_ac ratio
                analysis_TRAP(f_ac, voltage_times, r_mean)
            
            elif is_size_known == False and calibration_method == "average":
                f_TRAP, voltage_times, columns = prepare_analyzed_data_for_TRAP(force_values, voltage_map) # filter particles, produce time points array
                f_ac = f_TRAP / Phys_properties['calibr'] # calculate f_ac ratio
                r_TRAP_row = sizes.iloc[columns].to_numpy().transpose() / 2
                r_TRAP = np.tile(r_TRAP_row,(num_voltages,1))
                analysis_TRAP(forces_TRAP, voltage_times, r_TRAP)
    
            elif is_size_known == True and calibration_method == "heatmap":     
                f_ac = \
                    apply_heatmap_calibration(force_values, x_up_values, y_up_values)
                f_ac_TRAP, voltage_times, columns = prepare_analyzed_data_for_TRAP(f_ac, voltage_map) # filter particles, produce time points array            
                analysis_TRAP(f_ac_TRAP, voltage_times, r_mean)
                
            elif is_size_known == False and calibration_method == "heatmap":    
                f_ac = \
                    apply_heatmap_calibration(force_values, x_up_values, y_up_values)
                f_ac_TRAP, voltage_times, columns = prepare_analyzed_data_for_TRAP(f_ac, voltage_map) # filter particles, produce time points array                        
                r_TRAP_row = sizes.iloc[columns].to_numpy().transpose() / 2
                r_TRAP = np.tile(r_TRAP_row,(num_voltages,1))
                analysis_TRAP(f_ac_TRAP, voltage_times, r_TRAP)
 
#%% 
if calibration_method == "heatmap" and is_calibration == True:
    os.chdir(dir_calibr)
    produce_heatmap(combined_fxy_heatmap)
# QuantitativeAcoustophoresis
Python scripts to analyze acoustophoresis data (https://pubs.acs.org/doi/10.1021/acsnanoscienceau.2c00002)

To run the script, following Python libraries are needed:
math, matplotlib, nptdms, numpy, os, pandas, ptitprince, scipy 


Depending on the experimental conditions, certain parameters have to be changed. The example below showcases the values for a calibration measurement (average force) of 5.31 um PS beads:

temperature = 20                    # Celsius degrees; room temperature or the one set using a special stage
objective_magn = 20                 # 10/20/40 depending on which objective is used
diameter_sample = 5.31              # microbead size

is_size_known = True                # True - for monodisperse microbeads, False for other samples
is_temp_dep = False                 # True - use the temperature to calculate precisely the viscosity of the medium
is_density_known = False            # True - the sample density is known, False - otherwise
is_calibration = True               # True - for calibration sample, False - otherwise
is_time_resolved = False
calibration_method = "average"      # "heatmap" or "average"       
density_method = "full"             # "linear", "full" or "forces" fit
sample = "PS"                       # material of the sample (PS/Silica/PMMA) or its type (MEF/PMN/GUV/GUV gel)
medium = "PBS"                      # PBS/water/glucose/HEPES+


To analyze PMN cell data after that, certain adjustments need to be made

temperature = 20                    # Celsius degrees; room temperature or the one set using a special stage
objective_magn = 20                 # 10/20/40 depending on which objective is used
diameter_sample = 5.31              # microbead size

is_size_known = False               # True - for monodisperse microbeads, False for other samples
is_temp_dep = False                 # True - use the temperature to calculate precisely the viscosity of the medium
is_density_known = False            # True - the sample density is known, False - otherwise
is_calibration = False              # True - for calibration sample, False - otherwise
is_time_resolved = False
calibration_method = "average"      # "heatmap" or "average"       
density_method = "full"             # "linear", "full" or "forces" fit
sample = "PMN"                      # material of the sample (PS/Silica/PMMA) or its type (MEF/PMN/GUV/GUV gel)
medium = "HEPES+"                    # PBS/water/glucose/HEPES+


Also, the paths to the folders containing the calibration data (usually PS data measured on the same day), experimental data and script should be put in the according lines
dir_calibr = r"E:\MyFiles\Experimental data\Acoustophoresis\20200923 - PMNs and empty liposomes 1 of 3\AP_4.65um_SAV_coated_PS_beads_10-0p45V_20x_recalibration"
dir_data = r"E:\MyFiles\Experimental data\Acoustophoresis\20200923 - PMNs and empty liposomes 1 of 3\AP_PMNs_HEPES+_10-0p65V_20x\3"
dir_script = r"E:\MyFiles\Scripts\Acoustophoresis"

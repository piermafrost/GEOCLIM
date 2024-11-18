'''
Examples of use of the function generating GEOCLIM oceanic inputs

> select the block you want to execute with variable "CASE"
'''

from BC_generator import oce_output_2_GEOCLIM_BC


CASE = 'None'



# ------------------ ##
# Examples for IPSL: ##
# ------------------ ##

# IPSL CTRL (old way, without mask):
if CASE == 'IPSL old PI ctrl':
    oce_output_2_GEOCLIM_BC(co2=[284.7],
                            input_files=['piControl_SE_2750_2849_1Y_grid_T.nc'],
                            root='../../INPUT/IPSL/CTRL_CM5A2/',
                            grid_file='../../INPUT/IPSL/coordinates_ORCA_R2_n.nc',
                            latitude='nav_lat', z='deptht', temperature='thetao',
                            horiz_weight='cell_area', x_weight='e1t', y_weight='e2t', z_weight='e3t',
                            outdir='../../INPUT/COMBINE/IPSL_CTRL_old/',
                            check_plot=True)

# IPSL CTRL with oce mask, single epicontinental basin (REF)
if CASE == 'IPSL PI ref':
    oce_output_2_GEOCLIM_BC(co2=[284.7],
                            input_files=['piControl_SE_2750_2849_1Y_grid_T.nc'],
                            input_u_files=['piControl_SE_2750_2849_1Y_grid_U.nc'],
                            input_v_files=['piControl_SE_2750_2849_1Y_grid_V.nc'],
                            input_w_files=['piControl_SE_2750_2849_1Y_grid_W.nc'],
                            grid_file='../../INPUT/IPSL/coordinates_ORCA_R2_n.nc',
                            root='../../INPUT/IPSL/CTRL_CM5A2/',
                            mask_nc_file='oce_masks_IPSL_piCTRL.nc', mask_nc_var='mask', mask_values=(1,0),
                            basins_kw={'lat_split': (-65.,),
                                       'depth_split': (100., 900.),
                                       'polar_auto_identification': True},
                            latitude='nav_lat', z='deptht', temperature='thetao',
                            u='uocetr_eff', v='vocetr_eff', w=None, w_method='divergence',
                            horiz_weight='cell_area', x_weight='e1t', y_weight='e2t', z_weight='e3t',
                            flux_correction='RLS', special_wrap='ORCA-T',
                            outdir='../../INPUT/COMBINE/IPSL-PI-ref/',
                            export_mask_fields=True)#, check_plot=True)

# IPSL CTRL with multiple basins, split epicontinental box
if CASE == 'IPSL PI AtlPac':
    oce_output_2_GEOCLIM_BC(co2=[284.7],
                            input_files=['piControl_SE_2750_2849_1Y_grid_T.nc'],
                            input_u_files=['piControl_SE_2750_2849_1Y_grid_U.nc'],
                            input_v_files=['piControl_SE_2750_2849_1Y_grid_V.nc'],
                            input_w_files=['piControl_SE_2750_2849_1Y_grid_W.nc'],
                            grid_file='../../INPUT/IPSL/coordinates_ORCA_R2_n.nc',
                            root='../../INPUT/IPSL/CTRL_CM5A2/',
                            mask_nc_file='oce_masks_IPSL_piCTRL.nc', mask_nc_var='mask_ArcAtlPacSO', mask_values=(1,2,3,4),
                            basins_kw={'lat_split': (),
                                       'depth_split': (100., 900.),
                                       'split_epicont': True,
                                       'polar_mask_values': (1,4),
                                       'polar_auto_identification': False},
                            latitude='nav_lat', z='deptht', temperature='thetao',
                            u='uocetr_eff', v='vocetr_eff', w=None, w_method='divergence',
                            horiz_weight='cell_area', x_weight='e1t', y_weight='e2t', z_weight='e3t',
                            flux_correction='RLS', special_wrap='ORCA-T',
                            outdir='../../INPUT/COMBINE/IPSL-PI_AtlPac/',
                            export_mask_fields=True)#, check_plot=True)

# IPSL 90Ma, single climate, no input mask:
# NOTE: netCDF inputs are not stored on current GitHub repository
if CASE == 'IPSL 90Ma':
    oce_output_2_GEOCLIM_BC(co2=[2*284.7],
                            input_files=['CPL-90Ma-ORB7-2X_SE_7250_7349_1Y_grid_T.nc'],
                            input_u_files=['CPL-90Ma-ORB7-2X_SE_7250_7349_1Y_grid_U.nc'],
                            input_v_files=['CPL-90Ma-ORB7-2X_SE_7250_7349_1Y_grid_V.nc'],
                            input_w_files=['CPL-90Ma-ORB7-2X_SE_7250_7349_1Y_grid_W.nc'],
                            root='../../INPUT/IPSL/90Ma_Laugie/2X/Orb7/',
                            grid_file='../../INPUT/IPSL/90Ma_Laugie/grid/PALEORCA2.90MaCorrected_grid.nc',
                            basins_kw={'lat_split': (-65.,58.)},
                            latitude='nav_lat', z='deptht', temperature='thetao',
                            u='uocetr_eff', v='vocetr_eff', w=None, w_method='divergence',
                            horiz_weight='srft', x_weight='e1t', y_weight='e2t', z_weight='e3t',
                            outdir='test_IPSL_90Ma_Orb7a/',
                            special_wrap='ORCA-F', # NOTE: F-pivot is used for paleo grid, whereas T-pivot is used for the standard pre-industrial grid
                            export_mask_fields=True)


# ----------------- ##
# Example for FOAM: ##
# ----------------- ##

if CASE == 'FOAM':
    oce_output_2_GEOCLIM_BC(co2=[280, 1120],
                            root='../../INPUT/FOAM/PIctrl_cpld/',
                            input_files=['CTRL.AP.ocean.ctrlAC_1368W.nc', 'CTRL.AP.ocean.ctrlAC_1368W_1120ppm.nc'],
                            grid_file='../../INPUT/FOAM/grid_FOAM_128x128.nc',
                            basins_kw={'lat_split': (-63.,57.)},
                            temperature='TEMP', u='U', v='V', w_method='divergence', flux_correction='Positive Least Square',
                            latitude='lat', z='lev', horiz_weight='area', x_weight='cell_x_len', y_weight='cell_y_len', z_weight='thickness',
                            outdir='../../INPUT/COMBINE/FOAM-CTRL-nomask/',
                            export_mask_fields=True)



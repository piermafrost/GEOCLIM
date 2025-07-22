from BC_generator import oce_output_2_GEOCLIM_BC


WHICH = '360Ma-split'


if WHICH == 'CTRL':
    oce_output_2_GEOCLIM_BC(co2=[280.],
                            input_files=['PIctrl_cpld/0tc_1368.0W_OrbPD_ocean_280ppm.nc'],
                            grid_file='grid_FOAM_128x128.nc',
                            root='../../INPUT/FOAM/',
                            latitude='lat', z='lev', temperature='TEMP',
                            mask_nc_file='PIctrl_cpld/mask_oce_FOAM_PIctrl.nc', mask_nc_var='mask_arctic', mask_values=(1,0),
                            basins_kw={'lat_split': (-63.,), 'depth_split': (100., 1000.), 'polar_auto_identification': True},
                            u='U', v='V', w=None, w_method='divergence',
                            horiz_weight='area', x_weight='cell_x_len', y_weight='cell_y_len', z_weight='thickness',
                            flux_correction='RLS',
                            outdir='../../INPUT/COMBINE/FOAM-CTRL/',
                            export_mask_fields=True, check_plot=True)


#if WHICH == '360Ma':
#    oce_output_2_GEOCLIM_BC(co2=[140., 280., 560., 1120., 2240., 4480., 8960.],
#                            input_files=['360Ma_cpld/360eb_1327.33W_EccN_ocean_140ppm.nc',
#                                         '360Ma_cpld/360eb_1327.33W_EccN_ocean_280ppm.nc',
#                                         '360Ma_cpld/360eb_1327.33W_EccN_ocean_560ppm.nc',
#                                         '360Ma_cpld/360eb_1327.33W_EccN_ocean_1120ppm.nc',
#                                         '360Ma_cpld/360eb_1327.33W_EccN_ocean_2240ppm.nc',
#                                         '360Ma_cpld/360eb_1327.33W_EccN_ocean_4480ppm.nc',
#                                         '360Ma_cpld/360eb_1327.33W_EccN_ocean_8960ppm.nc'],
#                            grid_file='grid_FOAM_128x128.nc',
#                            root='../../INPUT/FOAM/',
#                            latitude='lat', z='lev', temperature='TEMP',
#                            mask_nc_file=None,
#                            basins_kw={'lat_split': (-50.,70.,), 'depth_split': (100., 1000.), 'polar_auto_identification': True,
#                                       'epicont_depth': 250}, # <= there are only 2 levels in the paleobathymetry reconstruction: ~100m and ~250m
#                            u='U', v='V', w=None, w_method='divergence',
#                            horiz_weight='area', x_weight='cell_x_len', y_weight='cell_y_len', z_weight='thickness',
#                            flux_correction='RLS',
#                            outdir='../../INPUT/COMBINE/FOAM-360Ma/',
#                            export_mask_fields=True)


if WHICH == '360Ma-split':
    oce_output_2_GEOCLIM_BC(co2=[140., 280., 560., 1120., 2240., 4480., 8960.],
                            input_files=['360Ma_cpld/360eb_1327.33W_EccN_ocean_140ppm.nc',
                                         '360Ma_cpld/360eb_1327.33W_EccN_ocean_280ppm.nc',
                                         '360Ma_cpld/360eb_1327.33W_EccN_ocean_560ppm.nc',
                                         '360Ma_cpld/360eb_1327.33W_EccN_ocean_1120ppm.nc',
                                         '360Ma_cpld/360eb_1327.33W_EccN_ocean_2240ppm.nc',
                                         '360Ma_cpld/360eb_1327.33W_EccN_ocean_4480ppm.nc',
                                         '360Ma_cpld/360eb_1327.33W_EccN_ocean_8960ppm.nc'],
                            grid_file='grid_FOAM_128x128.nc',
                            root='../../INPUT/FOAM/',
                            latitude='lat', z='lev', temperature='TEMP',
                            mask_nc_file='360Ma_cpld/360eb_OceBasMsk.nc',
                            mask_nc_var='mask_open', mask_values=(1,2,3,4),
                            epicmask_nc_var='mask_epic', epicmask_values=(1,2,3,4,5),
                            basins_kw={'lat_split': (), 'depth_split': (100., 1000.),
                                       'polar_auto_identification': False, 'polar_mask_values': (4,),
                                       'epicont_depth': 250}, # <= there are only 2 levels in the paleobathymetry reconstruction: ~100m and ~250m
                            u='U', v='V', w=None, w_method='divergence',
                            horiz_weight='area', x_weight='cell_x_len', y_weight='cell_y_len', z_weight='thickness',
                            flux_correction='RLS',
                            outdir='../../INPUT/COMBINE/FOAM-360Ma-split/',
                            export_mask_fields=True)


from basinmap_editor import  make_routing_map

make_routing_map(['../../../INPUT/IPSL/90Ma_Laugie/other/basinmap_90Ma-ORB7a-2X.nc',
                  '../../../INPUT/IPSL/90Ma_Laugie/2X/Orb1/CPL-90Ma-ORB1-2X_SE_6950_7049_1Y_sechiba.nc',
                  '../../../INPUT/IPSL/90Ma_Laugie/grid/PALEORCA2.90MaCorrected_grid.nc'],
                 lon='lon', lat='lat', lsm='Contfrac',
                 oce_basin_mask='../../../INPUT/COMBINE/IPSL-90Ma-HR/basins_masks.nc',
                 ocelon='nav_lon', ocelat='nav_lat',
                 watershed_map='basinmap', outlets_map='nbrivers',
                 output_file='GEOCLIM_basin_map_IPSL-90Ma-HR.nc',
                 invert_axis='y')


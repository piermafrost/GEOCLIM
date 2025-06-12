from basinmap_editor import  make_routing_map


WHICH = '3bas Arct-3'


# Single global ocean, split epicontinental basin
# -----------------------------------------------

if WHICH == 'split epicont' or WHICH == 'ALL':
    make_routing_map(['../../../INPUT/IPSL/90Ma_Laugie/other/basinmap_90Ma-ORB7a-2X.nc',
                    '../../../INPUT/IPSL/90Ma_Laugie/2X/Orb1/CPL-90Ma-ORB1-2X_SE_6950_7049_1Y_sechiba.nc',
                    '../../../INPUT/IPSL/90Ma_Laugie/grid/PALEORCA2.90MaCorrected_grid.nc'],
                    lon='lon', lat='lat', lsm='Contfrac',
                    oce_basin_mask='../../../INPUT/COMBINE/IPSL-90Ma-splitepic/basins_masks.nc',
                    ocelon='nav_lon', ocelat='nav_lat',
                    watershed_map='basinmap', outlets_map='nbrivers',
                    output_file='GEOCLIM_basin_map_IPSL-90Ma-splitepic.nc', output_var='geoclim_routing_map_globoce',
                    invert_axis='y')


# Single global ocean, split epicontinental basin, Arctic ocean split in 3 vertical levels (instead of 2)
# -------------------------------------------------------------------------------------------------------

if WHICH == 'split epicont Arct-3' or WHICH == 'ALL':
    make_routing_map(['../../../INPUT/IPSL/90Ma_Laugie/other/basinmap_90Ma-ORB7a-2X.nc',
                    '../../../INPUT/IPSL/90Ma_Laugie/2X/Orb1/CPL-90Ma-ORB1-2X_SE_6950_7049_1Y_sechiba.nc',
                    '../../../INPUT/IPSL/90Ma_Laugie/grid/PALEORCA2.90MaCorrected_grid.nc'],
                    lon='lon', lat='lat', lsm='Contfrac',
                    oce_basin_mask='../../../INPUT/COMBINE/IPSL-90Ma-splitepic-Arct3-2X/basins_masks.nc',
                    ocelon='nav_lon', ocelat='nav_lat',
                    watershed_map='basinmap', outlets_map='nbrivers',
                    output_file='GEOCLIM_basin_map_IPSL-90Ma-splitepic-Arct3.nc', output_var='geoclim_routing_map_globoce_Arct3',
                    invert_axis='y')


# Split oceans Arctic - Atlantic - Pacific - South Pole Pac - South Pole Atl
# --------------------------------------------------------------------------

if WHICH == 'Atl-Pac' or WHICH == 'ALL':
    make_routing_map(['../../../INPUT/IPSL/90Ma_Laugie/other/basinmap_90Ma-ORB7a-2X.nc',
                    '../../../INPUT/IPSL/90Ma_Laugie/2X/Orb1/CPL-90Ma-ORB1-2X_SE_6950_7049_1Y_sechiba.nc',
                    '../../../INPUT/IPSL/90Ma_Laugie/grid/PALEORCA2.90MaCorrected_grid.nc'],
                    lon='lon', lat='lat', lsm='Contfrac',
                    oce_basin_mask='../../../INPUT/COMBINE/IPSL-90Ma-AtlPac/basins_masks.nc',
                    ocelon='nav_lon', ocelat='nav_lat',
                    watershed_map='basinmap', outlets_map='nbrivers',
                    output_file='GEOCLIM_basin_map_IPSL-90Ma-AtlPac.nc', output_var='geoclim_routing_map_AtlPac',
                    invert_axis='y')


# Split oceans Arctic - Atlantic - Pacific - South Pole Pac - South Pole Atl, Arctic ocean split in 3 vertical levels (instead of 2)
# ----------------------------------------------------------------------------------------------------------------------------------

if WHICH == 'Atl-Pac Arct-3' or WHICH == 'ALL':
    make_routing_map(['../../../INPUT/IPSL/90Ma_Laugie/other/basinmap_90Ma-ORB7a-2X.nc',
                    '../../../INPUT/IPSL/90Ma_Laugie/2X/Orb1/CPL-90Ma-ORB1-2X_SE_6950_7049_1Y_sechiba.nc',
                    '../../../INPUT/IPSL/90Ma_Laugie/grid/PALEORCA2.90MaCorrected_grid.nc'],
                    lon='lon', lat='lat', lsm='Contfrac',
                    oce_basin_mask='../../../INPUT/COMBINE/IPSL-90Ma-AtlPac-Arct3-2X/basins_masks.nc',
                    ocelon='nav_lon', ocelat='nav_lat',
                    watershed_map='basinmap', outlets_map='nbrivers',
                    output_file='GEOCLIM_basin_map_IPSL-90Ma-AtlPac-Arct3.nc', output_var='geoclim_routing_map_AtlPac_Arct3',
                    invert_axis='y')


# Split oceans Arctic - N Atlantic - S Atlantic - Pacific - South Pole Pac - South Pole Atl
# -----------------------------------------------------------------------------------------

if WHICH == '3bas' or WHICH == 'ALL':
    make_routing_map(['../../../INPUT/IPSL/90Ma_Laugie/other/basinmap_90Ma-ORB7a-2X.nc',
                    '../../../INPUT/IPSL/90Ma_Laugie/2X/Orb1/CPL-90Ma-ORB1-2X_SE_6950_7049_1Y_sechiba.nc',
                    '../../../INPUT/IPSL/90Ma_Laugie/grid/PALEORCA2.90MaCorrected_grid.nc'],
                    lon='lon', lat='lat', lsm='Contfrac',
                    oce_basin_mask='../../../INPUT/COMBINE/IPSL-90Ma-3bas/basins_masks.nc',
                    ocelon='nav_lon', ocelat='nav_lat',
                    watershed_map='basinmap', outlets_map='nbrivers',
                    output_file='GEOCLIM_basin_map_IPSL-90Ma-3bas.nc', output_var='geoclim_routing_map_3bas',
                    invert_axis='y')


# Split oceans Arctic - N Atlantic - S Atlantic - Pacific - South Pole Pac - South Pole Atl, Arctic ocean split in 3 vertical levels (instead of 2)
# -------------------------------------------------------------------------------------------------------------------------------------------------

if WHICH == '3bas Arct-3' or WHICH == 'ALL':
    make_routing_map(['../../../INPUT/IPSL/90Ma_Laugie/other/basinmap_90Ma-ORB7a-2X.nc',
                    '../../../INPUT/IPSL/90Ma_Laugie/2X/Orb1/CPL-90Ma-ORB1-2X_SE_6950_7049_1Y_sechiba.nc',
                    '../../../INPUT/IPSL/90Ma_Laugie/grid/PALEORCA2.90MaCorrected_grid.nc'],
                    lon='lon', lat='lat', lsm='Contfrac',
                    oce_basin_mask='../../../INPUT/COMBINE/IPSL-90Ma-3bas-Arct3/basins_masks.nc',
                    ocelon='nav_lon', ocelat='nav_lat',
                    watershed_map='basinmap', outlets_map='nbrivers',
                    output_file='GEOCLIM_basin_map_IPSL-90Ma-3bas-Arct3.nc', output_var='geoclim_routing_map_3bas_Arct3',
                    invert_axis='y')



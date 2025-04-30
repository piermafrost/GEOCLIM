from basinmap_editor import make_routing_map

make_routing_map(input_files=('../../../INPUT/FOAM/360Ma_cpld/360eb_lsm.nc',
  '/home/piermafrost/Documents/Crisis-project_FOAM_Feifei/boundary_conditions/Topobathy_360eb_postslarti_cor.nc'),
                 lon='lon', lat='lat', lsm='lsm', watershed_map='TOPO',
                 oce_basin_mask='../../../INPUT/COMBINE/FOAM-360Ma-split/basins_masks.nc',
                 output_file='../../../INPUT/FOAM/360Ma_cpld/GEOCLIM_basin_map_FOAM360Ma-split.nc')


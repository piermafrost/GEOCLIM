"ARTIFICIAL" MODIFICATION OF INPUT DATA:

The original sedimenting areas of intermediate boxes were zero, because seafloor is 3500m deep everywhere on FOAM bathymetry
(except along the coasts, where it is 50m or 250m deep).

Because GEOCLIM cannot handle null sedimenting area, a small non-zero area has been set, by "transferring" area from deep boxes to intermediate boxes.
Try `diff oce_surf.dat oce_surf.ORIGINAL.dat` and `diff surf_sedi.dat surf_sedi.ORIGINAL.dat` to see the modifications.

"ARTIFICIAL" MODIFICATION OF INPUT DATA:

The original sedimenting areas of intermediate boxes were zero, because seafloor is 3500m deep everywhere on FOAM bathymetry
(except along the coasts, where it is 50m or 250m deep).
Because GEOCLIM cannot handle null sedimenting area, a small non-zero area has been set, by "transferring" area from deep boxes to intermediate boxes.
Try `diff oce_surf.dat oce_surf.ORIGINAL.dat` and `diff surf_sedi.dat surf_sedi.ORIGINAL.dat` to see the modifications.
Concerned boxes are #: 5, 6, 8, 9, 10, 11

A second consequences is that, in "box_sf_bnd_len.dat" costal deep boxes were not connected to open intermediated boxes, and intermediate boxes were not conected to open deep boxes.
Therefore, the "connexion" length "costal deep --> open deep" have been copied to "costal deep --> open intermediate" and to "open intermediate --> open deep"
Try `diff box_sf_bnd_len.dat box_sf_bnd_len.ORIGINAL.dat` to see the modifications.
Concerned boxes are #: 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21

lastly, northern high-latitude boxes had to be manually declared as "polar", and because there is an epicontinental box specific for southern high-latitude, it has been declared "polar" to reduce its primary productivity.
Try `diff indice_polar.dat indice_polar.ORIGINAL.dat` to see the modifications.
Concerned boxes are #: 1, 2, 3, 20, 21


# astylo


+ bio
	+ read_fits
	+ write_fits
	+ ext_WCS
	+ read_hdf5
	+ write_hdf5
	+ read_ascii
	+ read_csv
	+ write_csv

+ plot
	+ plotool
		+ figure
		+ set_border
		+ Cartesian2d
		+ set_ax
		+ plot
		+ set_font
		+ fill
		+ text
		+ save
		+ show
	+ plot2d
	+ plot2d_m

+ lib
	+ fclean
	+ pix2sr
	+ sr2arcsec2
	+ rad2arcsec
	+ celest2deg
	+ f_lin
	+ f_lin1
	+ f_lin0
	+ gaussian
	+ gaussian2D
	+ rms
	+ nanrms
	+ std
	+ closest
	+ bsplinterpol
	+ MCerror

+ <sup>__&dagger;__</sup>proc
	+ wclean
	+ interfill
	+ hextract
	+ improve
		+ addunc
		+ slice
		+ rectangle
	+ slicube
		+ image
		+ wave
		+ slice_names
	+ crop
		+ image
		+ wave
	+ iproject
		+ image
		+ wave
		+ slice_names
	+ iconvolve
		+ spitzer_irs
		+ choker
		+ do_conv
		+ image
		+ wave
		+ slice_names
	+ sextract
	+ imontage

+ calib
	+ intercalib
		+ synthetic_photometry
	+ spec2phot
	+ phot2phot
	+ specorrect
	+ photometry_profile

<sup>__&dagger;__</sup> This module calls modules from the same package.


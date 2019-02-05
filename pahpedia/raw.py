import numpy as np
import hdunit as hd


fits_image_path = "../data/m83/"

fits_image_filename = 'm83_ll2_cube'
fits_ref0_filename = 'm83_ll2_cube_ref'

#data = hd.rdata3d(fits_image_path, fits_image_filename)[0][0]
#hd.cphdr(data, fits_image_path, fits_image_filename, fits_image_path, fits_ref0_filename)


fits_raw_filename = 'm83_akari_cube_raw'
fits_image_filename = 'm83_akari_cube'
fits_unc_filename = 'm83_akari_cube_unc'
fits_ref0_filename = 'm83_akari_cube_ref'

data0 = hd.rakari(fits_image_path, fits_raw_filename)[0]
data1 = hd.rakari(fits_image_path, fits_raw_filename)[1]
hd.cd2pc(data0[100,:,:], fits_image_path, fits_raw_filename, fits_image_path, fits_ref0_filename)
hd.cd2pc(data0, fits_image_path, fits_raw_filename, fits_image_path, fits_image_filename)
hd.cd2pc(data1, fits_image_path, fits_raw_filename, fits_image_path, fits_unc_filename)

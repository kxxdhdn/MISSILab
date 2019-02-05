; Read the file in and assign it to the 'param' variable;
; assign the header information to variables for use later.
param = READ_CSV('kern.csv', N_TABLE_HEADER=1)
;HELP, param, /STRUCTURES
im_file = param.FIELD1
kern_file = param.FIELD2
;print, im_file[0]
N = N_elements(im_file)

do_we_write = 0
FOR I=0,N-1 DO BEGIN
	im = MRDFITS("../data/m83/slices/" + im_file[I] + ".fits",0,hdr)
	kern = MRDFITS("../data/kernels/" + kern_file[I] + ".fits",0,hdrkern)
	do_the_convolution,im,hdr,kern,hdrkern,imconv,hdrconv,$
    	               result_kernel_image,result_kernel_header,do_we_write
	MWRFITS, imconv, "../data/m83/convolved/"+ im_file[I] + "_conv.fits", hdrconv, /CREATE
ENDFOR

END


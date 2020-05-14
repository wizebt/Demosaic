function bl(I,ColorID)
	"""
	Demosaic color patterns in 3 channels
	I is gray scale image from bayer mosaic
	rgb is a 3 D array with color channel
	"""
	#if ColorID<8 || ColorID>11 return end
	h,w=size(I)
	rgb=zeros(Int,3,h,w)
	r1,c1,r2,c2=1:2:h,1:2:w,2:2:h,2:2:w
	#Sensor patterns RGGB, GRBG, GBRG, BGGR
	sp=[1 2 2 3;2 1 3 2;2 3 1 2;3 2 2 1]
	i=ColorID-7
	rgb[sp[i,1],r1,c1]=I[r1,c1]
	rgb[sp[i,1],r1,2:2:w-1]=(I[r1,1:2:w-2]+I[r1,3:2:w])/2
	rgb[sp[i,1],2:2:h-1,c1]=(I[1:2:h-2,c1]+I[3:2:end,c1])/2
	rgb[:,:,1]
	#rgb[sp[i,2],r1,c2]=I[r1,c2]
	#rgb[sp[i,3],r2,c1]=I[r2,c1]
	#rgb[sp[i,4],r2,c2]=I[r2,c2]

	return rgb
end

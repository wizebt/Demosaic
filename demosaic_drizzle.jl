"""
Drizzle demosaic for BAYER patterns split RGB channels
(C) Emmanuel Brandt 2019-20
(8, RGGB),(9,GRBG),(10,GBRG),(11,BGGR)
"""

function demosaic_drizzle(I, ColorID)

	rgb = zeros(eltype(I), size(I,1), size(I,2), 3)
	idx = [1 2 2 3;2 1 3 2;2 3 1 2;3 2 2 1]
	c = ColorId - 7
	
	rgb[1:2:end,1:2:end,idx[c,1]] = I[1:2:end,1:2:end] 
	rgb[1:2:end,2:2:end,idx[c,2]] = I[1:2:end,2:2:end] 
	rgb[2:2:end,1:2:end,idx[c,3]] = I[2:2:end,1:2:end] 
	rgb[2:2:end,2:2:end,idx[c,4]] = I[2:2:end,2:2:end] 

	return rgb
end

function demosaic_drizzle1(I)
	h, w = size(I)
	rgb=zeros(eltype(I),h,w,3)
	
	rgb[1:2:h, 1:2:w, 1] = I[1:2:h, 1:2:w]
	rgb[1:2:h, 2:2:w, 2] = I[1:2:h, 2:2:w]
	rgb[2:2:h, 1:2:w, 2] = I[2:2:h, 1:2:w]
	rgb[2:2:h, 2:2:w, 3] = I[2:2:h, 2:2:w]

	rgb
end

function demosaic_drizzle2(I)
	h, w = size(I)
	rgb = zeros(eltype(I), h, w, 3)
	x, y = 1:2:h, 1:2:w
	rgb[x, y, 1] = I[x, y]
	x1, y1 = 2:2:h, 2:2:w
	rgb[x, y1, 2] = I[x, y1]
	rgb[x1, y, 2] = I[x1, y]
	rgb[x1, y1, 3] = I[x1, y1]

	rgb
end

"""
	median9bayer(I::Array{Float32,2})

Median 3x3 filter for bayer image.
"""
function median55bayer(I::Array{Float32,2})
	return mapwindow(median9, I, (3,3))
end


"""
	median55bayer(I::Array{Float32,2})

Median 3x3 filter for bayer image.
"""
function median55bayer(I::Array{Float32,2})
	return mapwindow(median9, I, (3,3))
end

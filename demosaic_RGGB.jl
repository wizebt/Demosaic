"""
	demosaic_RGGB(I::Array{Float32,2})

Fast bilinear demosaic of CFA image with RGGB bayer sensor .
"""
function demosaic_RGGB(I::Array{Float32,2})
    h, w = size(I)
    rgb = zeros(Float32, h, w, 3)
    a, b = 0.5f0, 0.25f0

    @inbounds for j = 1:2:w-1, i = 1:2:h-1
        i1, j1 = i + 1, j + 1
        rgb[i, j, 1] = I[i, j]
        rgb[i, j1, 2] = I[i, j1]
        rgb[i1, j, 2] = I[i1, j]
        rgb[i1, j1, 3] = I[i1, j1]
    end
    @inbounds for j = 1:2:w-1, i = 2:2:h-2
        rgb[i, j, 1] = (I[i-1, j] + I[i+1, j])a
    end
    @inbounds for j = 2:2:w-2, i = 1:2:h-1
        rgb[i, j, 1] = (I[i, j-1] + I[i, j+1])a
    end
    @inbounds for j = 2:2:w-1, i = 2:2:h-1
        rgb[i, j, 1] = (I[i-1, j-1] + I[i-1, j+1] + I[i+1, j-1] + I[i+1, j-1])b
        rgb[i, j, 2] = (I[i-1, j] + I[i, j-1] + I[i, j+1] + I[i+1, j])b
    end
    @inbounds for j = 3:2:w-1, i = 3:2:h-1
        rgb[i, j, 2] = (I[i-1, j] + I[i, j-1] + I[i, j+1] + I[i+1, j])b
    end
    @inbounds for j = 2:2:w, i = 3:2:h-1
        rgb[i, j, 3] = (I[i-1, j] + I[i+1, j])a
    end
    @inbounds for j = 3:2:w-1, i = 2:2:h
        rgb[i, j, 3] = (I[i, j-1] + I[i, j+1])a
    end
    @inbounds for j = 3:2:w, i = 3:2:h-1
        rgb[i, j, 3] = (I[i-1, j-1] + I[i-1, j+1] + I[i+1, j-1] + I[i+1, j+1])b
    end

    return rgb
end

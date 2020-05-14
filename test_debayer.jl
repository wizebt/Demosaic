using ImageView, Images, OffsetArrays, ComputationalResources

function t(; debayer::String="drizzle")
	I = rand(Float32, 128, 256)
	
	if debayer == "drizzle"
		println("Drizzle")
		rgb = dd(I,8)
	elseif debayer == "bilinear"
		println("Bilinear")
		rgb = demosaic_bilinear(I,8)
	else
		error("error $debayer")
	end
	
	myshow(rgb)
	return nothing
end	

using ImageView
function myshow(I)
	imshow(colorview(RGB, permutedims(I, (3,1,2))))
end

"""
	dd(I::Array{Float32,2}, ColorID::Int) -> rgb

Drizzle Debayerdevectorized.
"""
function dd(I::Array{Float32,2}, ColorID)
    h, w = size(I)
    rgb = zeros(Float32, size(I, 1), size(I, 2), 3)
    idx = [1 2 2 3; 2 1 3 2; 2 3 1 2; 3 2 2 1]

    if 7 < ColorID < 12
        cid = ColorID - 7
    else
        error("LuckyWizard>dd : ColorID $ColorID not supported")
    end
    a, b, c, d = idx[:, cid]

    @inbounds for j = 1:2:w, i = 1:2:h
        i1, j1 = i + 1, j + 1
        
        rgb[i, j, a] = I[i, j]
        rgb[i, j1, b] = I[i, j1]
        rgb[i1, j, c] = I[i1, j]
        rgb[i1, j1, d] = I[i1, j1]
    end

    return rgb
end

using ImageFiltering
"""
	debayer(I, ColorID)
        
Demosaic bilinear color pattern. 
I CFA image
RGB is a 3 D array with color channel
ColorID is the sensor pattern code
"""
function demosaic_bilinear(I::Array{Float32,2}, ColorID)
    if ColorID == 8 # RGGB
        ind = [
            1 1    # Red
            2 1    # Green 1
            1 2    # Green 2
            2 2
        ] # Blue
    elseif ColorID == 9 # GBRG
        ind = [
            1 2    # Red
            1 1    # Green 1
            2 2    # Green 2
            2 1
        ]  # Blue
    elseif ColorID == 10 # GBRG
        ind = [
            2 1    # Red
            1 1    # Green 1
            2 2    # Green 2
            1 2
        ]  # Blue
    elseif ColorID == 11 # BAYER_BGGR
        ind = [
            2 2    # Red
            2 1    # Green 1
            1 2    # Green 2
            1 1
        ]  # Blue
    else
        error("demosaic: Invalid ColorID ", ColorID)
    end

    # Image processing kernels
    krb = [ 0.25 0.5 0.25
            0.5 1 0.5
            0.25 0.5 0.25 ]
    kg = [  0 0.25 0
            0.25 1 0.25
            0 0.25 0  ]

    # Pad input image by 1 pixel to avoid artifacts at image borders
    I = padarray(I, Pad(:reflect, 1, 1))
    z = 0.0f0
    # Shifts starting indices of R/G/B channels
    ind = -ind .+ 3 #ind=3*ones(Int16,4,2)-indrgb
    h, w = size(I)
    I = OffsetArray(I, 1:h, 1:w)
    rgb = zeros(Float32, h, w, 3)
    y1, x1 = 2:h-1, 2:w-1

    # Red Channel
    channel = zeros(Float32, h, w)
    y, x = ind[1, 1]:2:h, ind[1, 2]:2:w
    channel[y, x] = I[y, x]
    channel = imfilter(CPU1(Algorithm.FIR()), channel, krb)
    channel[y, x] = I[y, x]
    rgb[y1, x1, 1] = channel[y1, x1]

    # Green Channel (both pixels)
    fill!(channel, z)
    y, x = ind[2, 1]:2:h, ind[2, 2]:2:w
    channel[y, x] = I[y, x]
    y, x = ind[3, 1]:2:h, ind[3, 2]:2:w
    channel[y, x] = I[y, x]
    channel = imfilter(CPU1(Algorithm.FIR()), channel, kg)
    rgb[y1, x1, 2] = channel[y1, x1]

    # Blue Channel
    fill!(channel, z)
    y, x = ind[4, 1]:2:h, ind[4, 2]:2:w
    channel[y, x] = I[y, x]
    channel = imfilter(CPU1(Algorithm.FIR()), channel, krb)
    rgb[y1, x1, 3] = channel[y1, x1]

    return rgb[y1, x1, :]
end

"""
	demosaic_bilinear(I::Array{Float32,2}, ColorID) -> rgb
        
Demosaic bilinear color pattern. 
- `I`: CFA image
- `ColorID`: Sensor pattern code
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
    krb = centered(	[0.25 0.5 0.25; 0.5 1 0.5; 0.25 0.5 0.25 ])
    kg = centered(	[0 0.25 0; 0.25 1 0.25; 0 0.25 0 ])

    # Pad input image by 1 pixel to avoid artifacts at image borders
    #I = padarray(I, Pad(:reflect, 1, 1))
    z = 0.0f0
    # Shifts starting indices of R/G/B channels
    ind = -ind .+ 3 # ind=3*ones(Int16,4,2)-indrgb
    h, w = size(I)
    I = OffsetArray(I, 1:h, 1:w)
    #rgb = zeros(Float32, h, w, 3)
    rgb = Array{Float32,3}(undef, h, w, 3)
    y1, x1 = 2:h-1, 2:w-1

    # Red Channel
    channel = zeros(Float32, h, w)
    y, x = ind[1, 1]:2:h, ind[1, 2]:2:w
    channel[y, x] = I[y, x]
    #channel = imfilter(CPU1(Algorithm.FIR()), channel, krb)
    channel = imfilter(CPUThreads(Algorithm.FIR()), channel, krb)
    #channel = imfilter(CUDALibs(Algorithm.FIR()), channel, krb)
    #channel = imfilter(ArrayFireLibs(Algorithm.FIR()), channel, krb)
    #channel = imfilter(channel, krb)
    channel[y, x] = I[y, x]
    rgb[y1, x1, 1] = channel[y1, x1]

    # Green Channel (both pixels)
    fill!(channel, z)
    y, x = ind[2, 1]:2:h, ind[2, 2]:2:w
    channel[y, x] = I[y, x]
    y, x = ind[3, 1]:2:h, ind[3, 2]:2:w
    channel[y, x] = I[y, x]
    #channel = imfilter(CPU1(Algorithm.FIR()), channel, kg)
    channel = imfilter(CPUThreads(Algorithm.FIR()), channel, kg)
    #channel = imfilter(CUDALibs(Algorithm.FIR()), channel, kg)
    #channel = imfilter(ArrayFireLibs(Algorithm.FIR()), channel, kg)
    #channel = imfilter(channel, kg)
    rgb[y1, x1, 2] = channel[y1, x1]

    # Blue Channel
    fill!(channel, z)
    y, x = ind[4, 1]:2:h, ind[4, 2]:2:w
    #channel[y, x] = I[y, x]
    #channel = imfilter(CPU1(Algorithm.FIR()), channel, krb)
    channel = imfilter(CPUThreads(Algorithm.FIR()), channel, krb)
    #channel = imfilter(CUDALibs(Algorithm.FIR()), channel, krb)
    #channel = imfilter(ArrayFireLibs(Algorithm.FIR()), channel, krb)
    #channel = imfilter(channel, krb)
    rgb[y1, x1, 3] = channel[y1, x1]

    return rgb[y1, x1, :]
end

using Images, OffsetArrays
I=reshape(Float32.(1:24),4,6);

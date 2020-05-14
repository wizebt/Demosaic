using ImageFiltering, OffsetArrays

function debayer(I, ColorID)
    """
    Demosaic color patterns in 3 channels
    I CFA image
    RGB is a 3 D array with color channel
    ColorID is the sensor pattern code
    """

    #sensorAlignment = tolower(sensorAlignment);
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
        println("demosaic: Invalid ColorID ", ColorID)
    end

    ## Image processing kernels
    krb = [
        0.25 0.5 0.25
        0.5 1 0.5
        0.25 0.5 0.25
    ]
    kg = [
        0 0.25 0
        0.25 1 0.25
        0 0.25 0
    ]

    # Pad input image by 1 pixel to avoid artifacts at image borders
    I = padarray(I, Pad(:reflect, 1, 1))
    ## Shifts starting indices of R/G/B channels
    ind = -ind .+ 3 #ind=3*ones(Int16,4,2)-indrgb
    h, w = size(I)
    I = OffsetArray(I, 1:h, 1:w)
    rgb = zeros(h, w, 3)
    y1, x1 = 2:h-1, 2:w-1

    ## Red Channel
    channel = zeros(h, w)
    y, x = ind[1, 1]:2:h, ind[1, 2]:2:w
    channel[y, x] = I[y, x]
    channel = imfilter(channel, krb)
    rgb[y1, x1, 1] = channel[y1, x1]

    # Green Channel (both pixels)
    channel = zeros(h, w)
    y, x = ind[2, 1]:2:h, ind[2, 2]:2:w
    channel[y, x] = I[y, x]
    y, x = ind[3, 1]:2:h, ind[3, 2]:2:w
    channel[y, x] = I[y, x]
    channel = imfilter(channel, kg)
    rgb[y1, x1, 2] = channel[y1, x1]

    # Blue Channel
    channel = zeros(h, w)
    y, x = ind[4, 1]:2:h, ind[4, 2]:2:w
    channel[y, x] = I[y, x]
    channel = imfilter(channel, krb)
    rgb[y1, x1, 3] = channel[y1, x1]

    #return map(x->round(UInt16,x), rgb[y1,x1,:])
    rgb[y1, x1, :]
end

function demosaic_rggb(I::Array{Float32,2})
    h, w = size(I)
    rgb = zeros(Float32, h, w, 3)

    @inbounds for j = 1:2:w-1, i = 1:2:h-1
        i1, j1 = i + 1, j + 1
        rgb[i, j, 1] = I[i, j]
        rgb[i, j1, 2] = I[i, j1]
        rgb[i1, j, 2] = I[i1, j]
        rgb[i1, j1, 3] = I[i1, j1]
    end
    @inbounds for j = 1:2:w-1, i = 2:2:h-2
        rgb[i, j, 1] = 0.5 * (I[i-1, j] + I[i+1, j])
    end
    @inbounds for j = 2:2:w-2, i = 1:2:h-1
        rgb[i, j, 1] = 0.5 * (I[i, j-1] + I[i, j+1])
    end
    @inbounds for j = 2:2:w-1, i = 2:2:h-1
        rgb[i, j, 1] =
            0.25 * (I[i-1, j-1] + I[i-1, j+1] + I[i+1, j-1] + I[i+1, j-1])
        rgb[i, j, 2] = 0.25 * (I[i-1, j] + I[i, j-1] + I[i, j+1] + I[i+1, j])
    end
    @inbounds for j = 3:2:w-1, i = 3:2:h-1
        rgb[i, j, 2] = 0.25 * (I[i-1, j] + I[i, j-1] + I[i, j+1] + I[i+1, j])
    end
    @inbounds for j = 2:2:w, i = 3:2:h-1
        rgb[i, j, 3] = 0.5 * (I[i-1, j] + I[i+1, j])
    end
    @inbounds for j = 3:2:w-1, i = 2:2:h
        rgb[i, j, 3] = 0.5 * (I[i, j-1] + I[i, j+1])
    end
    @inbounds for j = 3:2:w, i = 3:2:h-1
        rgb[i, j, 3] =
            0.25 * (I[i-1, j-1] + I[i-1, j+1] + I[i+1, j-1] + I[i+1, j+1])
    end

    return rgb
end

# Test data
I = reshape(Float32.(1:48), 6, 8)

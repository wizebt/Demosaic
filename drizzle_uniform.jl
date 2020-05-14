function drizzle_uniform(I)
    h, w = size(I)
    rgb = Array{Float32,3}(undef, h, w, 3)

    for j = 1:2:w, i = 1:2:h
        rgb[i, j, 1] =
            rgb[i, j+1, 1] = rgb[i+1, j, 1] = rgb[i+1, j+1, 1] = I[i, j]
        rgb[i, j, 2] =
            rgb[i, j+1, 2] =
                rgb[i+1, j, 2] =
                    rgb[i+1, j+1, 2] = 0.5 * (I[i, j+1] + I[i+1, j])
        rgb[i, j, 3] =
            rgb[i, j+1, 3] = rgb[i+1, j, 3] = rgb[i+1, j+1, 3] = I[i+1, j+1]
    end

    return rgb
end

 function M = rescale_binary_map (im, scale)

    M = zeros(ceil(size(im,1)/scale), ceil(size(im,2)/scale));

    for i = scale:size(im,1)
        for j = scale: size(im,2)
            if(im(i,j)==1)
                M(int32(i/scale), int32(j/scale)) = 1;
            end
        end
    end

    M = logical(M);
 end
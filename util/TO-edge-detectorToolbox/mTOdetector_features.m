function f = mTOdetector_features(img, sigma, n, islogist)
% now only one channel: lumination

if (size(img,3)>1) %convert from rgb to gray
    img = rgb2gray(img);
end
img = double(img);

% compute gaussian derivation response map for each scale
for s = 1:length(sigma)
    Ix(:,:,s)   = resampled_filter_2d(img, @Gx_2d_op, sigma(s), n);
    Iy(:,:,s)   = resampled_filter_2d(img, @Gy_2d_op, sigma(s), n);
    grad_mag(:,:,s) = sqrt(Ix(:,:,s).^2+Iy(:,:,s).^2);
end

[d1, d2, d3] = size(grad_mag);
grad_mag = reshape(grad_mag, [d1*d2, d3]);


f = grad_mag';
if(islogist)
    f= [ones(1, d1*d2); f];
end


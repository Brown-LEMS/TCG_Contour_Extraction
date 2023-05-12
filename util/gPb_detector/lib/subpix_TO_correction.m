function [ TO_edgemap, TO_orientation, edginfo ] = subpix_TO_correction( edgemap,theta,thresh,sigma)
% Use NMS to detect the edgels - this NMS routine ouputs subpixel edgel tokens
margin = 3;
dirx = cos(mod(theta, pi));
diry = sin(mod(theta, pi));
[subpix_x, subpix_y, subpix_dir_x, subpix_dir_y, subpix_grad_x, subpix_grad_y]...
    = NMS_token(dirx, diry,edgemap, edgemap>thresh, margin);

% magnitude of the gradient at the maxima
mag_e = sqrt(subpix_grad_x.^2+subpix_grad_y.^2);

Hx = edgemap.*dirx;
Hy = edgemap.*diry;
% sigma = 2;
Px   = resampled_filter_2d(edgemap, @Gx_2d_op, sigma, 0);
Py   = resampled_filter_2d(edgemap, @Gy_2d_op, sigma, 0);
Pxx  = resampled_filter_2d(edgemap, @Gxx_2d_op, sigma, 0);
Pxy  = resampled_filter_2d(edgemap, @Gxy_2d_op, sigma, 0);
Pyy  = resampled_filter_2d(edgemap, @Gyy_2d_op, sigma, 0);

Hxx  = resampled_filter_2d(Hx, @Gx_2d_op, sigma, 0);
Hyy  = resampled_filter_2d(Hy, @Gy_2d_op, sigma, 0);
Hxy  = resampled_filter_2d(Hy, @Gx_2d_op, sigma, 0);

[yd, xd] = size(Hx);
[xx,yy] = meshgrid(1:xd,1:yd);
Px_e   = interp2(xx, yy, Px,  subpix_x, subpix_y, 'cubic');
Py_e   = interp2(xx, yy, Py,  subpix_x, subpix_y, 'cubic');
Pxx_e  = interp2(xx, yy, Pxx,  subpix_x, subpix_y, 'cubic');
Pyy_e  = interp2(xx, yy, Pyy,  subpix_x, subpix_y, 'cubic');
Pxy_e  = interp2(xx, yy, Pxy,  subpix_x, subpix_y, 'cubic');

Hx_e  = interp2(xx, yy, Hx, subpix_x, subpix_y, 'cubic');
Hy_e  = interp2(xx, yy, Hy, subpix_x, subpix_y, 'cubic');
Hxx_e = interp2(xx, yy, Hxx, subpix_x, subpix_y, 'cubic');
Hxy_e = interp2(xx, yy, Hxy, subpix_x, subpix_y, 'cubic');
Hyy_e = interp2(xx, yy, Hyy, subpix_x, subpix_y, 'cubic');

% compute [Fx, Fy] at zero crossings (second order equations on P = |hist-grad I|)
Fx_e = Px_e.*Hxx_e + Pxx_e.*Hx_e + Py_e.*Hxy_e + Pxy_e.*Hy_e;
Fy_e = Px_e.*Hxy_e + Pxy_e.*Hx_e + Py_e.*Hyy_e + Pyy_e.*Hy_e;

%normalize
F_mag = sqrt(Fx_e.^2 + Fy_e.^2);
Fx_e = Fx_e./F_mag;
Fy_e = Fy_e./F_mag;

% Edge Direction (tangent to the level set) is orthogonal to the gradient
subpix_dir_x2 = (-Fy_e);
subpix_dir_y2 = (Fx_e);

% construct the edge map
edginfo = [subpix_x' subpix_y' atan2(subpix_dir_y2', subpix_dir_x2') mag_e']; 
TO_edgemap = zeros(size(edgemap));
TO_orientation = zeros(size(edgemap));
for q = 1:size(edginfo, 1)
    x = max(round(edginfo(q, 1)),1); x =min(x,size(edgemap,2));
    y = max(round(edginfo(q, 2)),1); y =min(y,size(edgemap,1));
    TO_edgemap(y,x) = edginfo(q, 4);
    TO_orientation(y,x) = edginfo(q, 3);
end
        
% debug only
% figure;
% disp_edg(edginfo, 'b');
% axis([100 180 100 180]);

end


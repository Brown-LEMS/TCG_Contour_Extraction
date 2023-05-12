function TO_edges = subpix_TO_correction_gPb( gPb_orient,im,thresh,r)

norient=8;
gtheta = [1.5708    1.1781    0.7854    0.3927   0    2.7489    2.3562    1.9635];

[~, maxo] = max(gPb_orient,[],3);
gPb_ori = gtheta(maxo); %orientation of maxima

% nonmax suppression and max over orientations
[h,w,unused] = size(im);
gPb = zeros(h,w);
theta = zeros(h,w);
% r = 2.5;
% r = sigma;

for i = 1:norient,
  mask = (maxo == i);
  a = fitparab_old(gPb_orient(:,:,i),r,r,gtheta(i));
  gPbi = nonmax(max(0,a),gtheta(i));
  gPb = max(gPb,gPbi.*mask);
  theta = theta.*~mask + gtheta(i).*mask;
end
gPb_max = max(0,min(1,gPb));


%[unused, maxo] = max(gPb_orient,[],3);
ori_max = gtheta(maxo); %orientation of maxima

% three observations a,b,c around the max to fit the parabola
a = zeros(size(maxo));
b = zeros(size(maxo));
c = zeros(size(maxo));
for i=1:norient,
    a = a + gPb_orient(:,:,i).*(maxo==i);
    b = b + gPb_orient(:,:,i).*(mod(maxo-1, norient)==i);
    c = c + gPb_orient(:,:,i).*(mod(maxo+1, norient)==i);
end

wedgesize = pi/norient;

d = (b + c - 2 * a);
degen = (abs(d) < 1e-3); % degenerate mask
x = (wedgesize/2)*(b - c)./d;
gPb = ~degen.*(a + x.*(c - b)/(2*wedgesize) + x.*x.*d/(2*wedgesize*wedgesize)) + degen.*a; % if degenerate, just use the discrete one
theta = ~degen.*(mod(ori_max + x + 2*pi, pi)) + degen.*ori_max;

%gPb_orient = max(gPb_orient, [], 3);
gPb = max(0,min(1,gPb));
gPb_orient = gPb;

% Use NMS to detect the edgels - this NMS routine ouputs subpixel edgel tokens
margin = 3;
dirx = cos(mod(theta+pi/2, pi));
diry = sin(mod(theta+pi/2, pi));
[subpix_x, subpix_y, subpix_dir_x, subpix_dir_y, subpix_grad_x, subpix_grad_y] = NMS_token(dirx, diry,gPb_orient, gPb_orient>thresh, margin);

% gPb_orient = gPb_max;
% margin = 3;
% dirx = cos(mod(ori_max+pi/2, pi));
% diry = sin(mod(ori_max+pi/2, pi));
% [subpix_x, subpix_y, subpix_dir_x, subpix_dir_y, subpix_grad_x, subpix_grad_y] = NMS_token(dirx, diry,gPb_orient, gPb_orient>thresh, margin);

% magnitude of the gradient at the maxima
mag_e = sqrt(subpix_grad_x.^2+subpix_grad_y.^2);

Hx = gPb_orient.*dirx;
Hy = gPb_orient.*diry;
sigma = 0.7;
Px   = resampled_filter_2d(gPb_orient, @Gx_2d_op, sigma, 0);
Py   = resampled_filter_2d(gPb_orient, @Gy_2d_op, sigma, 0);
Pxx  = resampled_filter_2d(gPb_orient, @Gxx_2d_op, sigma, 0);
Pxy  = resampled_filter_2d(gPb_orient, @Gxy_2d_op, sigma, 0);
Pyy  = resampled_filter_2d(gPb_orient, @Gyy_2d_op, sigma, 0);

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
gPb_interp = interp2(xx, yy, gPb_orient, subpix_x, subpix_y, 'cubic');

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
%edg  = [subpix_x' subpix_y' atan2(subpix_dir_y2', subpix_dir_x2') mag_e'];

% construct the edge map
TO_edges  = [subpix_x' subpix_y' atan2(subpix_dir_y2', subpix_dir_x2') gPb_interp'];


%%%%%%%%%%% gPb

ind=find(gPb_max);
row=mod(ind,size(gPb_max,1));
row(row==0)=size(gPb_max,1);
sub=ind-row;
col=sub/size(im,1);
col=col+ones(size(col));
%col(row==50)=col(row==50)-1;
%gPb_max = max(gPb_orient,[],3);

mat=[];
for i=1:size(row,1)
mat=[mat; col(i) row(i) ori_max(row(i),col(i)) gPb_max(row(i),col(i))];
end

%mag_e = sqrt(mat(:,1).^2+mat(:,2).^2);
gPb_edges = mat;

end
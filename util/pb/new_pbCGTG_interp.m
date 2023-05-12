function [pb,theta] = new_pbCGTG_interp(im,beta,radius,norient)
% function [pb,theta] = pbCGTG(im,radius,norient)
% 
% Compute probability of boundary using CG and TG.
%
% David R. Martin <dmartin@eecs.berkeley.edu>
% April 2003

if nargin<3, radius=[0.01 0.02 0.02 0.02]; end
if nargin<4, norient=8; end
if numel(radius)==1, radius=radius*ones(1,4); end

% % beta from logistic fits (trainCGTG.m)
% if all(radius==[0.01 0.02 0.02 0.02]), % 64 textons
%   beta = [   -2.1593726e+00   1.3343505e+00   2.0588519e-02  -9.5044157e-02   2.7115263e-02y];
%   fstd = [    1.0000000e+00   4.0475308e-01   1.3302566e-01   2.0786224e-01   1.8890904e-01];
%   beta = beta ./ fstd;
% else
%   error(sprintf('no parameters for radius=[%g %g]\n',radius(1),radius(2)));
% end

% get gradients
[cg,tg,gtheta] = detCGTG(im,radius,norient);



% compute oriented pb
[h,w,unused] = size(im);
pball = zeros(h,w,norient);
for i = 1:norient,
  l = cg(:,:,1,i); l = l(:);
  a = cg(:,:,2,i); a = a(:);
  b = cg(:,:,3,i); b = b(:);
  t = tg(:,:,i); t = t(:);
  x = [ones(size(b)) l a b t];
  pbi = 1 ./ (1 + (exp(-x*beta')));
  pball(:,:,i) = reshape(pbi,[h w]);
end

% nonmax suppression and max over orientations
[unused,maxo] = max(pball,[],3);
pb = zeros(h,w);
theta = zeros(h,w);
r = 2.5;
for i = 1:norient,
  mask = (maxo == i);
  a = fitparab(pball(:,:,i),r,r,gtheta(i));
  pbi = nonmax(max(0,a),gtheta(i));
  pb = max(pb,pbi.*mask);
  theta = theta.*~mask + gtheta(i).*mask;
end
pb = max(0,min(1,pb));

%interperate the orientation
new_theta = zeros(h,w);
for i = 1:h
   for j = 1:w
%        if(pb(i,j)>0)  
           if(maxo(i,j)==1 || maxo(i,j)==norient)
               new_theta(i,j) = gtheta(maxo(i,j));
           else
                x = (maxo(i,j)-2:maxo(i,j))/norient*pi;
                y = [pball(i,j,maxo(i,j)-1), pball(i,j,maxo(i,j)),pball(i,j,maxo(i,j)+1)];
                s0 = 3; s1 = sum(x); s2 = sum(x.^2); s3 = sum(x.^3); s4 = sum(x.^4);
                A = [s4,s3,s2;s3,s2,s1;s2,s1,s0];
                d = [sum(x.^2.*y);sum(x.*y);sum(y)];
                parA = A\d;
                new_theta(i,j) = -parA(2)/2/parA(1);
           end
%        end
   end
end
theta = new_theta;

% mask out 1-pixel border where nonmax suppression fails
pb(1,:) = 0;
pb(end,:) = 0;
pb(:,1) = 0;
pb(:,end) = 0;





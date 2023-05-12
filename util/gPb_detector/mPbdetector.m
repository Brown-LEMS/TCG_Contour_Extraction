function [f] = mPbdetector(im)
dim = size(im,3);

%% mPb
if dim == 1,
   im(:,:,2)=im(:,:,1);im(:,:,3)=im(:,:,1); 
end
% get gradients
tic;
[bg1, bg2, bg3, cga1, cga2, cga3, cgb1, cgb2, cgb3, tg1, tg2, tg3, textons] = det_mPb(im);
fprintf('Local cues: %g\n', toc);

% smooth cues
gtheta = [1.5708    1.1781    0.7854    0.3927   0    2.7489    2.3562    1.9635];
tic;
filters = make_filters([3 5 10 20], gtheta);
for o = 1 : size(tg1, 3),
    bg1(:,:,o) = fitparab(bg1(:,:,o),3,3/4,gtheta(o),filters{1,o});
    bg2(:,:,o) = fitparab(bg2(:,:,o),5,5/4,gtheta(o),filters{2,o});
    bg3(:,:,o) = fitparab(bg3(:,:,o),10,10/4,gtheta(o),filters{3,o});

    cga1(:,:,o) = fitparab(cga1(:,:,o),5,5/4,gtheta(o),filters{2,o});
    cga2(:,:,o) = fitparab(cga2(:,:,o),10,10/4,gtheta(o),filters{3,o});
    cga3(:,:,o) = fitparab(cga3(:,:,o),20,20/4,gtheta(o),filters{4,o});

    cgb1(:,:,o) = fitparab(cgb1(:,:,o),5,5/4,gtheta(o),filters{2,o});
    cgb2(:,:,o) = fitparab(cgb2(:,:,o),10,10/4,gtheta(o),filters{3,o});
    cgb3(:,:,o) = fitparab(cgb3(:,:,o),20,20/4,gtheta(o),filters{4,o});

    tg1(:,:,o) = fitparab(tg1(:,:,o),5,5/4,gtheta(o),filters{2,o});
    tg2(:,:,o) = fitparab(tg2(:,:,o),10,10/4,gtheta(o),filters{3,o});
    tg3(:,:,o) = fitparab(tg3(:,:,o),20,20/4,gtheta(o),filters{4,o});

end
fprintf('Cues smoothing: %g\n', toc);


bg1 = max(bg1, [], 3);
bg2 = max(bg2, [], 3);
bg3 = max(bg3, [], 3);
cga1 = max(cga1, [], 3);
cga2 = max(cga2, [], 3);
cga3 = max(cga3, [], 3);
cgb1 = max(cgb1, [], 3);
cgb2 = max(cgb2, [], 3);
cgb3 = max(cgb3, [], 3);
tg1 = max(tg1, [], 3);
tg2 = max(tg2, [], 3);
tg3 = max(tg3, [], 3);

bg1 = bg1(:);
bg2 = bg2(:);
bg3 = bg3(:);
cga1 = cga1(:);
cga2 = cga2(:);
cga3 = cga3(:);
cgb1 = cgb1(:);
cgb2 = cgb2(:);
cgb3 = cgb3(:);
tg1 = tg1(:);
tg2 = tg2(:);
tg3 = tg3(:);

% if(dim==3)
    f = [ones(size(bg1)), bg1, bg2, bg3, cga1, cga2, cga3, cgb1, cgb2, cgb3, tg1, tg2, tg3]';

% [cg,tg] = detCGTG(im);
% l = max(squeeze(cg(:,:,1,:)),[],3);
% a = max(squeeze(cg(:,:,2,:)),[],3);
% b = max(squeeze(cg(:,:,3,:)),[],3);
% t = max(tg,[],3);
% l = l(:);
% a = a(:);
% b = b(:);
% t = t(:);
% f = [ ones(size(b)) l a b t ]';

%%
function filters = make_filters(radii, gtheta)

d = 2; 

filters = cell(numel(radii), numel(gtheta));
for r = 1:numel(radii),
    for t = 1:numel(gtheta),
        
        ra = radii(r);
        rb = ra / 4;
        theta = gtheta(t);
        
        ra = max(1.5, ra);
        rb = max(1.5, rb);
        ira2 = 1 / ra^2;
        irb2 = 1 / rb^2;
        wr = floor(max(ra, rb));
        wd = 2*wr+1;
        sint = sin(theta);
        cost = cos(theta);
        
        % 1. compute linear filters for coefficients
        % (a) compute inverse of least-squares problem matrix
        filt = zeros(wd,wd,d+1);
        xx = zeros(2*d+1,1);
        for u = -wr:wr,
            for v = -wr:wr,
                ai = -u*sint + v*cost; % distance along major axis
                bi = u*cost + v*sint; % distance along minor axis
                if ai*ai*ira2 + bi*bi*irb2 > 1, continue; end % outside support
                xx = xx + cumprod([1;ai+zeros(2*d,1)]);
            end
        end
        A = zeros(d+1,d+1);
        for i = 1:d+1,
            A(:,i) = xx(i:i+d);
        end
        
        % (b) solve least-squares problem for delta function at each pixel
        for u = -wr:wr,
            for v = -wr:wr,
                ai = -u*sint + v*cost; % distance along major axis
                bi = u*cost + v*sint; % distance along minor axis
                if (ai*ai*ira2 + bi*bi*irb2) > 1, continue; end % outside support
                yy = cumprod([1;ai+zeros(d,1)]);
                filt(v+wr+1,u+wr+1,:) = A\yy;
            end
        end
        
        filters{r,t}=filt;
    end
end

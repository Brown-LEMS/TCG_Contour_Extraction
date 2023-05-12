function bry_map_trans = transform_bry_map_given_GT_depth(ID1,ID2, bry_map1, datasetFilePath)

addpath(datasetFilePath);
ImgID1 = ID1;
ImgID2 = ID2;

% I1 = (single(rgb2gray(imread([datasetFilePath 'Camera', num2str(ImgID1) ,'.png']))) ./ 255);
% I2 = (single(rgb2gray(imread([datasetFilePath, 'Camera', num2str(ImgID2) ,'.png']))) ./ 255);
[h, w, ~] = size(bry_map1);
%II = [I1,255 .* ones(500,5),I2,255 .* ones(500,5),I3];
%figure;
%imagesc(II);
%axis off
%colormap gray
%Load Depth Map
depth1 = importdata([datasetFilePath 'Camera', num2str(ImgID1) ,'depth.txt']);
depth1 = transpose(reshape(depth1,500,500));
%depth1(find(depth1 >= 8)) = 8;
depth1 = flip(depth1);

% depth2 = importdata([datasetFilePath 'Camera',num2str(ImgID2) ,'depth.txt']);
% depth2 = transpose(reshape(depth2,500,500));
% %depth2(find(depth2 >= 8)) = 8;
% depth2 = flip(depth2);

%Load Camera Parameters
camPara = importdata([datasetFilePath 'Camera',num2str(ImgID1),'Param.txt']);
K1 = camPara(1:9);
K1 = reshape(K1,[3,3]);
K1 = K1';
P1 = camPara(10:21);
P1 = reshape(P1,[4,3]);
P1 = P1';
RT1 = inv(K1) * P1;
R1 = RT1(1:3,1:3); 
T1 = RT1(1:3,4); 

% Rw1 = R1'; 
% Cw1 = -Rw1 * T1;
% P1 = K1 * [Rw1' -Rw1' * Cw1];

camPara2 = importdata([datasetFilePath 'Camera',num2str(ImgID2),'Param.txt']);
K2 = camPara2(1:9);
K2 = reshape(K2,[3,3]);
P2 = camPara2(10:21);
K2 = K2';
P2 = reshape(P2,[4,3]);
P2 = P2';
RT2 = inv(K2) * P2;
R2 = RT2(1:3,1:3);
T2 = RT2(1:3,4);
% Rw2 = R2';
% Cw2 = -Rw2 * T2;

R12  = R2 * inv(R1);
T12 = - (R2 * inv(R1)) * T1 + T2;

% %% Check points in 3D
% CorresImgPixel1 = zeros(4, h*w);
% x_coords = repmat(1:w, h, 1);
% y_coords = repmat((1:h)', 1, w);
% 
% CorresImgPixel1(1, :) = x_coords(:);
% CorresImgPixel1(2, :) = y_coords(:);
% CorresImgPixel1(3, :) = I1(sub2ind([h,w], y_coords(:), x_coords(:)));
% CorresImgPixel1(4, :) = depth1(sub2ind([h,w], y_coords(:), x_coords(:)));


%% Remap Depth Map
x_coords = repmat(1:w, h, 1);
y_coords = repmat((1:h)', 1, w);

p1 = [x_coords(:) y_coords(:)];


p1homo = inv(K1) * [p1'; ones(1, size(p1,1))];
% z = depth1(sub2ind([h,w], round(p1(:,2))+1,round(p1(:,1)+1)));
% look for the foreground depth in a local nbr as the depth of an edge
z_vec = [];
for dx = [-1 0 1]
    x_coords = round(p1(:,1))+1 + dx;
    x_coords =  max(x_coords, 1);
    x_coords = min(x_coords, w);
    for dy = [-1 0 1]
        y_coords = round(p1(:,2))+1 + dy;
        y_coords =  max(y_coords, 1);
        y_coords = min(y_coords, h);

        z_vec = [z_vec depth1(sub2ind([h,w], y_coords, x_coords))];
    end
end
z = min(z_vec, [], 2);
rho1 = z'./sqrt( (p1homo(1,:).^2 + p1homo(2,:).^2 + 1));

%% transform edges

p23D = R12 * (repmat(rho1, 3,1) .* p1homo) + repmat(T12, 1, size(p1homo, 2)) ;
p2 = K2 * p23D;
p2 = p2./ repmat(p2(3, :), 3, 1);
p2(2,:) = min(round(p2(2,:)), h);
p2(2,:) = max(round(p2(2,:)), 1);
p2(1,:) = min(round(p2(1,:)), w);
p2(1,:) = max(round(p2(1,:)), 1);

bry_map_trans = zeros(size(bry_map1));
bry_map_trans(sub2ind([h,w], p2(2,:), p2(1,:))) = bry_map1(:);



end
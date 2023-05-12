function [edges_trans, cfrags_trans, cfrags_idx_trans] = transform_Contours_given_GT_depth_v2(ID1,ID2, cfrags_idx1, edges1, datasetFilePath)
%%%%%%%%%%%%%%%%%%% version 2 %%%%%%%%%%%%%%%%%
%%%%%%%% a. remove the connection to junction if it is transformed to
%%%%%%%% a surface with huge gap in depth
%%%%%%%% b. remove the potion of contours that is occluded after
%%%%%%%% transformation

addpath(datasetFilePath);
ImgID1 = ID1;
ImgID2 = ID2;

I1 = (single(rgb2gray(imread([datasetFilePath 'Camera', num2str(ImgID1) ,'.png']))) ./ 255);
I2 = (single(rgb2gray(imread([datasetFilePath, 'Camera', num2str(ImgID2) ,'.png']))) ./ 255);
[h, w, ~] = size(I1);
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

depth2 = importdata([datasetFilePath 'Camera',num2str(ImgID2) ,'depth.txt']);
depth2 = transpose(reshape(depth2,500,500));
%depth2(find(depth2 >= 8)) = 8;
depth2 = flip(depth2);

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

p1 = edges1(:, 1:2);
p1homo = inv(K1) * [p1'; ones(1, size(edges1,1))];
% z = depth1(sub2ind([h,w], round(p1(:,2))+1,round(p1(:,1)+1)));
% look for the foreground depth in a local nbr as the depth of an edge
z_vec = [];
for dx = -1:1
    x_coords = round(p1(:,1))+1 + dx;
    x_coords =  max(x_coords, 1);
    x_coords = min(x_coords, w);
    for dy = -1:1
        y_coords = round(p1(:,2))+1 + dy;
        y_coords =  max(y_coords, 1);
        y_coords = min(y_coords, h);

        z_vec = [z_vec depth1(sub2ind([h,w], y_coords, x_coords))];
        
        if(dx ==0 && dy==0)
            z_m = depth1(sub2ind([h,w], y_coords, x_coords));
        end
    end
end
z1 = min(z_vec, [], 2);
% z_ind = find(abs(z_m-min(z_vec, [], 2))<0.1);
% z1(z_ind) = z_m(z_ind);
rho1 = z1'./sqrt( (p1homo(1,:).^2 + p1homo(2,:).^2 + 1));

%% transform edges

p23D = R12 * (repmat(rho1, 3,1) .* p1homo) + repmat(T12, 1, size(p1homo, 2)) ;

p2 = K2 * p23D;
p2 = p2./ repmat(p2(3, :), 3, 1);

rho2 =  p23D(3,:);
p2homo = p23D./ repmat(p23D(3, :), 3, 1);
z2 = rho2 .* sqrt( p2homo(1,:).^2 + p2homo(2,:).^2 +1);


edges_trans = edges1;
edges_trans(:, 1:2) = p2(1:2,:)';
edges_trans(:, 1) =  max(edges_trans(:, 1), 0);
edges_trans(:, 1) =  min(edges_trans(:, 1), w-1);
edges_trans(:, 2) =  max(edges_trans(:, 2), 0);
edges_trans(:, 2) =  min(edges_trans(:, 2), h-1);

%%Get Contours
cfrags_trans = cell(1,length(cfrags_idx1));
cfrags_idx_trans =  cfrags_idx1;
for i = 1:length(cfrags_idx1)
   cur_edges = edges_trans(cfrags_idx1{i}, :);
   
   %%%%%%% check artifact junctions, and remove them
   if(length(cfrags_idx1{i})<=1)
       cfrags_trans{i} = [];
       cfrags_idx_trans{i} = [];
       continue;
   end
   if( sqrt(sum((cur_edges(1, 1:2) - cur_edges(2, 1:2)).^2)) > 5 )
%    if( abs(z2(cfrags_idx1{i}(1)) - z2(cfrags_idx1{i}(2))) > 5 )
       cur_edges(1,:) = [];
       cfrags_idx1{i}(1) = [];
   end
   if(length(cfrags_idx1{i})<=1)
       cfrags_trans{i} = [];
       cfrags_idx_trans{i} = [];
       continue;
   end
   if( sqrt(sum((cur_edges(end, 1:2) - cur_edges(end-1, 1:2)).^2)) > 5 )
%    if( abs(z2(cfrags_idx1{i}(end)) - z2(cfrags_idx1{i}(end-1))) > 5 )
       cur_edges(end,:) = [];
       cfrags_idx1{i}(end) = [];
   end
   if(length(cfrags_idx1{i})<=1)
       cfrags_trans{i} = [];
       cfrags_idx_trans{i} = [];
       continue;
   end

   %%%%%%% locate occluded portion of a contour and break the contour
   z_cur = z2(cfrags_idx1{i});
   x_cur = round(cur_edges(:,1))+1;
   x_cur = max(x_cur, 1);
   x_cur = min(x_cur, w);
   y_cur = round(cur_edges(:,2))+1;
   y_cur = max(y_cur, 1);
   y_cur = min(y_cur, w);
   depth_2_gt = depth2(sub2ind([h,w], y_cur, x_cur));
   
%    is_occlude = (z_cur' - depth_2_gt)>5;
   occlude_idx = find((z_cur' - depth_2_gt)>0.5);
   if(~isempty(occlude_idx))
        cur_edges(occlude_idx, :) = [];
        cfrags_idx1{i}(occlude_idx) = [];
       
%        keyboard;
   end
   
   cfrags_trans{i} = cur_edges;
   cfrags_idx_trans{i} = cfrags_idx1{i};
end

cfrags_trans = cfrags_trans(~cellfun('isempty',cfrags_trans));
cfrags_idx_trans = cfrags_idx_trans(~cellfun('isempty',cfrags_idx_trans)); 

end
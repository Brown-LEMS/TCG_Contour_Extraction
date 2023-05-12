function [correspContours] = GetGroundTruthCorrespondingCurves(ID1,ID2,CEM1, datasetFilePath)

addpath(datasetFilePath);
ImgID1 = ID1;
ImgID2 = ID2;

I1 = (single(rgb2gray(imread([datasetFilePath 'Camera', num2str(ImgID1) ,'.png']))) ./ 255);
I2 = (single(rgb2gray(imread([datasetFilePath, 'Camera', num2str(ImgID2) ,'.png']))) ./ 255);
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

%% Check points in 3D
counter = 1;
for i = 1:size(I1,1)
    for j = 1:size(I1,2)
        CorresImgPixel1(1,counter) = j;
        CorresImgPixel1(2,counter) = i;
        CorresImgPixel1(3,counter) = I1(i,j);
        CorresImgPixel1(4,counter) = depth1(i,j);
        
%         CorresImgPixel2(1,counter) = j;
%         CorresImgPixel2(2,counter) = i;
%         CorresImgPixel2(3,counter) = I2(i,j);
%         CorresImgPixel2(4,counter) = depth2(i,j);
 
        counter = counter + 1;
    end
end
%% Remap Depth Map
for i = 1:size(CorresImgPixel1,2)
    p1 = CorresImgPixel1(1:2,i);
    % p1home is the 3D location of p1 in camera-cented coordinates
    p1homo = inv(K1) * [p1;1];
    % depth1 record the distance from carmara center C1 to the 3-D point
    % need to convert to the distance along the Z axis of the camera
    % centered coodinates: rho1 = depth1 / sqrt(1 + p1homo.x^2 = p1homo.y^2)
    z = depth1(p1(2),p1(1));
    rho1(i) = sqrt(z^2 ./ (p1homo(1).^2 + p1homo(2).^2 + 1));
    rhoMap1(p1(2),p1(1)) = rho1(i);
end

%%Get Contours
contoursInOneView = CEM1{2};


%% Get Correspondence
% R12 = inv(Rw2) * Rw1;
% T12 = inv(Rw2) * (Cw1-Cw2);
for i = 1:length(contoursInOneView)
    %Interpolate Depth
    contourPoints = contoursInOneView{i}(:,1:2);
    [X,Y] = meshgrid(1:500);
    interpedDepth = interp2(X,Y,rhoMap1,contourPoints(:,1) + 1,contourPoints(:,2) + 1);
    correspContour = [];
    for j = 1:size(contourPoints,1)
        p1 = transpose(contourPoints(j,1:2));
        d1 = interpedDepth(j);
        
        p1 = [p1+ 1; 1] ;
        p1homo = inv(K1) * p1;
        % when Z in p1homo is always 1, d1*p1homo scale it at the right
        % 3D location in camera1-centered coordinates
        % rotation by R12, follow by transition T12, transform it to the 3D
        % locaiton in camera2-centered coordinates, 
        p23D = R12 * (d1 * p1homo) + T12 ;
%         % normalize by Z to cast it to the depth of image scale 
%         p2homo = p23D ./ p23D(3);
%         % transform to location in camera 2 image
%         p2 = K2 * p2homo;
        
        p2 = K2 * p23D;
        p2 = p2./p2(3);
        
        correspContour(j,:) = [p2(1) - 1,p2(2) - 1];
    end
    correspContours{i} = correspContour;
end
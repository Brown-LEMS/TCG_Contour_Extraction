inPath  = '../../data/singleImage';
inFolder = inPath;
outPath = '../../results/singleImage';
if ~exist(outPath,'dir')
    mkdir(outPath);
end
imgList= dir([inPath '/*.PNG']);
numImage = length(imgList);
for iI= 1:numImage
    tic
    im = imread([inPath '/' imgList(iI).name]);
    if (size(im,3)==3)
        im= rgb2gray(im);
    end
    img_name = imgList(iI).name(1:end-4);
    
    % Find edges using the Canny operator with default
    % edgeim = edge(im,'canny',[],5);
	edgeim = edge(im,'canny');

    figure(1), imshow(1-edgeim);
    %set(gcf,'PaperPositionMode','auto');
    %filename = [outPath '/' imgList(iI).name(1:end-4)  '_canny.png'];
    %saveas(gcf,filename);
    
    % Link edge pixels together into lists of sequential edge points, one
    % list for each edge contour. A contour/edgelist starts/stops at an 
    % ending or a junction with another contour/edgelist.
    % Here we discard contours less than 18 pixels long.

    [edgelist, labelededgeim] = edgelink(edgeim, 19);
    
    % Display the edgelists with random colours for each distinct edge 
    % in figure 2
    h=drawedgelist(edgelist, size(im), 2, 'rand', 5); axis off 

    filename = [outPath '/' imgList(iI).name(1:end-4)  '_edgelink.png'];
    saveas(gcf,filename);
    toc
    close all
end

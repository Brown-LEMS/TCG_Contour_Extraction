function disp_edg(edg, col)
%function disp_edg(edg, col)
%
% Display the edgemap on the current figure in the given color
%   The edgels are displayed as little line segments.
%
% (c) LEMS, Brown University
% Amir Tamrakar (amir_tamrakar@brown.edu)
% October 2007

hold on;

[m,n] = size(edg);

% plot the edgels
for i=1:m,
    plot([edg(i,1)+0.5*cos(edg(i,3)) edg(i,1)-0.5*cos(edg(i,3))], [edg(i,2)+0.5*sin(edg(i,3)) edg(i,2)-0.5*sin(edg(i,3))], col);
    %plot(edg(i,1), edg(i,2), 'r.');
end

hold off;
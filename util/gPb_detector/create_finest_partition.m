function [ws_wt] = create_finest_partition(pb_oriented)

pb = max(pb_oriented,[],3);
ws = watershed(pb);
ws_bw = (ws == 0);

contours = fit_contour(double(ws_bw));
angles = zeros(numel(contours.edge_x_coords), 1);

for e = 1 : numel(contours.edge_x_coords)
    if contours.is_completion(e), continue; end
    v1 = contours.vertices(contours.edges(e, 1), :);
    v2 = contours.vertices(contours.edges(e, 2), :);

    if v1(2) == v2(2),
        ang = 90;
    else
        ang = atan((v1(1)-v2(1)) / (v1(2)-v2(2)));
    end
    angles(e) = ang*180/pi;
end

orient = zeros(numel(contours.edge_x_coords), 1);
orient((angles<-78.75) | (angles>=78.75)) = 1;
orient((angles<78.75) & (angles>=56.25)) = 2;
orient((angles<56.25) & (angles>=33.75)) = 3;
orient((angles<33.75) & (angles>=11.25)) = 4;
orient((angles<11.25) & (angles>=-11.25)) =5;
orient((angles<-11.25) & (angles>=-33.75)) = 6;
orient((angles<-33.75) & (angles>=-56.25)) = 7;
orient((angles<-56.25) & (angles>=-78.75)) = 8;

ws_wt = zeros(size(ws_bw));
for e = 1 : numel(contours.edge_x_coords)
    if contours.is_completion(e), continue; end
    for p = 1 : numel(contours.edge_x_coords{e}),
        ws_wt(contours.edge_x_coords{e}(p), contours.edge_y_coords{e}(p)) = ...
            max(pb_oriented(contours.edge_x_coords{e}(p), contours.edge_y_coords{e}(p), orient(e)), ws_wt(contours.edge_x_coords{e}(p), contours.edge_y_coords{e}(p)));
    end
    v1=contours.vertices(contours.edges(e,1),:);
    v2=contours.vertices(contours.edges(e,2),:);
    ws_wt(v1(1),v1(2))=max( pb_oriented(v1(1),v1(2), orient(e)),ws_wt(v1(1),v1(2)));
    ws_wt(v2(1),v2(2))=max( pb_oriented(v2(1),v2(2), orient(e)),ws_wt(v2(1),v2(2)));
end
ws_wt=double(ws_wt);
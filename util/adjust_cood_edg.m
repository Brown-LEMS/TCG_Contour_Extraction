function adjust_cood_edg (edg_file_in, edg_file_out, adjust)
    addpath io/

    [edges, edgemap] = load_edg(edg_file_in);

    edges(:,1) = edges(:,1) + adjust;
    edges(:,2) = edges(:,2) + adjust;

    [h,w] = size(edgemap);

    save_edg(edg_file_out, edges, [w,h]);

end

function [cfrags, cfrags_idx] = prune_noise_curves(cfrags, cfrags_idx, edgemap_soft, params)

% remove extremely short and low avg contrast curves
% remove duplicate paths between two graph nodes
% GOAL: remove those short curves affecting merging
len_th = params.noise_len_th;
prob_th = params.noise_prob_th;
[h,w] = size(edgemap_soft);
G_test = construct_fac_graph_from_curve_fragments (cfrags_idx, cfrags);

cid_to_remove = [];
for vid =1:length(G_test.var)
   
    % consider to prune free-end curves
    if(length(G_test.var(vid).nbrs_fac)==1)


        cid = G_test.var(vid).nbrs_fac;
        cur_c = cfrags{cid};
        c_len = contour_length_mex(cur_c');    

%         probs = edgemap_soft(sub2ind([h,w], round(cur_c(:,2))+1, round(cur_c(:,1))+1));
%         avg_prob = mean(probs);
        if(c_len<len_th) %&& avg_prob<prob_th)
            cid_to_remove = [cid_to_remove cid];
%         else if(c_len<2*len_th && avg_prob<prob_th)
%             cid_to_remove = [cid_to_remove cid];                
        end
    % look for the case to short connected curves connect to two ends of a third short curve and make a cycle    
    elseif (length(G_test.var(vid).nbrs_fac)==2)
        cid1 = G_test.var(vid).nbrs_fac(1);
        cid2 = G_test.var(vid).nbrs_fac(2);
        if(cid1 == cid2) % skip a circle end
            continue;
        end
        if(length(cfrags_idx{cid1})>len_th || length(cfrags_idx{cid2})>len_th)
            continue;
        end
        
        % find child vids
        cvid1 = G_test.fac(cid1).nbrs_var(G_test.fac(cid1).nbrs_var~=vid);
        cvid2 = G_test.fac(cid2).nbrs_var(G_test.fac(cid2).nbrs_var~=vid);
        
        c_inter = intersect(G_test.var(cvid1).nbrs_fac, G_test.var(cvid2).nbrs_fac);
        if(~isempty(c_inter))
            cid_to_remove = [cid_to_remove cid1 cid2];
        end
        
    % look for duplicate links between two nodes    
    elseif(length(G_test.var(vid).nbrs_fac)==3)
        
        for j = 1:length(G_test.var(vid).nbrs_fac)
           cid = G_test.var(vid).nbrs_fac(j);
           if(length(cfrags_idx{cid})>2*len_th)
               continue;
           end
            
           cvid = G_test.fac(cid).nbrs_var(G_test.fac(cid).nbrs_var~=vid);
           if(isempty(cvid)) % skip the cycle
               continue;
           end
           
           for k = 1:length(G_test.var(cvid).nbrs_fac)
               ccid = G_test.var(cvid).nbrs_fac(k);
               if(ccid==cid)
                   continue;
               end
               
              if(length(cfrags_idx{ccid})>2*len_th)
                   continue;
              end
              
              ccvid = G_test.fac(ccid).nbrs_var(G_test.fac(ccid).nbrs_var~=cvid);
              
              % found the duplicate link
              if(ccvid==vid)
                  cur_c1 = cfrags{cid};
                  probs1 = edgemap_soft(sub2ind([h,w], round(cur_c1(:,2))+1, round(cur_c1(:,1))+1));

                  cur_c2 = cfrags{ccid};
                  probs2 = edgemap_soft(sub2ind([h,w], round(cur_c2(:,2))+1, round(cur_c2(:,1))+1));
                  
                  if(mean(probs1) < mean(probs2))
                      cid_to_remove = [cid_to_remove cid];
                  else
                      cid_to_remove = [cid_to_remove ccid];
                  end
                  
              end
           end 
           
        end
        
        
    end
    
end
% cid_to_remove = unique(cid_to_remove);
% cfrags(cid_to_remove) = [];
% cfrags_idx(cid_to_remove) = [];
% G_test = construct_fac_graph_from_curve_fragments (cfrags_idx, cfrags);
% cid_to_remove = [];
% remove individual noise curves
for cid =1:length(G_test.fac)
   
    cur_c_idx = cfrags_idx{cid};
    % remove tiny cycles
    if(cur_c_idx(1)==cur_c_idx(end))
        cur_c = cfrags{cid};
        c_len =contour_length_mex(cur_c');    
        if(c_len< 3*len_th)
            cid_to_remove = [cid_to_remove cid];            
        end
        continue;
    end
       
    % isolated curves
    nbr_vids = G_test.fac(cid).nbrs_var;
    if(length(G_test.var(nbr_vids(1)).nbrs_fac)==1 && length(G_test.var(nbr_vids(2)).nbrs_fac)==1)
        cur_c = cfrags{cid};
        c_len =contour_length_mex(cur_c');    
        probs = edgemap_soft(sub2ind([h,w], round(cur_c(:,2))+1, round(cur_c(:,1))+1));
        avg_prob = mean(probs);
        
        if(avg_prob<prob_th && c_len< 3*len_th)
            cid_to_remove = [cid_to_remove cid];    
        elseif(c_len < len_th) % delete extremely short ones directly
            cid_to_remove = [cid_to_remove cid];        
        end
        continue;
    end
    
        
end
cid_to_remove = unique(cid_to_remove);
cfrags(cid_to_remove) = [];
cfrags_idx(cid_to_remove) = [];

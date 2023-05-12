function cfrags_prob = compute_cfrags_prob_given_probmap (cfrags, ProbMap)

%%%%%%%% given cfrags is written in c++ coordinates
[h,w] = size(ProbMap);
cfrags_prob = zeros(1, length(cfrags));
for c = 1:length(cfrags)
   cur_frag = cfrags{c};
   if(size(cur_frag, 1)<2)
       cfrags_prob(c) = 0;
       continue;
   end
   x_coords = round(cur_frag(2:end-1,1))+1;
   y_coords = round(cur_frag(2:end-1,2))+1;
   x_coords = max(x_coords, 1);
   x_coords = min(x_coords, w);
   y_coords = max(y_coords, 1);
   y_coords = min(y_coords, h);   
   
   probs = ProbMap(sub2ind([h,w], y_coords, x_coords));
   cfrags_prob(c) = mean(probs);
end

end
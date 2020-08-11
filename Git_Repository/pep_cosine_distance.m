% This program calculates the cosine distance of calculated ratios from the
% (assumed) ratios ("known_ratios") and returns a column with cosine distances

function cosine_distance=pep_cosine_distance(measured_ratios,known_ratios)
rows=size(measured_ratios,1);
%initialize storage for peptide correlation
cosine_distance = zeros(rows,1);
% create peptide matrixes for each protein
for ndx=1:rows %loop through all peptides
      %calculate cosine distance for all peptides 
      if sum(measured_ratios(ndx,:)) > 0    %Check if ratio was calculated otherwise return 1 for distance
        cosine_distance(ndx)=(pdist([measured_ratios(ndx,:);known_ratios],'cosine'));    
      else
        cosine_distance(ndx)=1; 
      end
end
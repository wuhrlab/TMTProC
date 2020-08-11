%Convolve the isotopic precursor distribution k-1 times with the isotopic
%envelope of TMT alone, Then convolve resulting vector with Isotopic
%impurity matrix I(TMT)


function P_TMT=convolve_precursor_w_TMT(isotope_prec_wo_TMT,TMT_iso_imp_matrix,num_TMT)
%Define the TMT isotopic envelope vector t(TMT) for each of the used
%channels by summing up rows of the respective impurity matrices I(TMT)
TMT_iso_vector=sum(TMT_iso_imp_matrix,1);
temp_isotopic_envelope = isotope_prec_wo_TMT;
%Convolve num_TMT-1 times the precursor isotopic envelope with TMT
for index=1:(num_TMT-1)
    temp_isotopic_envelope = conv(temp_isotopic_envelope,TMT_iso_vector);
end
%Convolve the vector of the isotopic precursor minus 1TMT with the last TMT
%which results in the matrix P(TMT) defining position and fragmentation
%pattern
P_TMT = zeros(11,size(temp_isotopic_envelope,2)+2);
for index = 1:11
    %Calculate the Precursor_TMT matrix P(TMT)
    P_TMT(index,:) = conv(TMT_iso_imp_matrix(index,:), temp_isotopic_envelope);
end
%Only save columns -1 to +9 for Precursor TMT matrix P (12 columns)
P_TMT = P_TMT(:,num_TMT:num_TMT+11);
end
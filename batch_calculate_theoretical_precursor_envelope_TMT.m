function theo_precursor_w_TMT=batch_calculate_theoretical_precursor_envelope_TMT(exp)

theo_precursor_w_TMT=zeros(size(exp.data.z,1),8);
for index=1:size(exp.data.z,1)
if sum(exp.ratios(index,:)) > 0.5 %Check if ratios were calculated
    matrix_precursor=calc_TMT_position_and_isotope_distribution(exp.ratios(index,:),exp.data.theoretical_precursor_envelope(index,:),exp.data.num_TMT(index));
    vector_precursor=sum(matrix_precursor,1);
    theo_precursor_w_TMT(index,:)=vector_precursor;
end
end
end
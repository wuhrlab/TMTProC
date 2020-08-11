function num_TMT = calculate_num_TMT(pep_sequences)
length_array = size(pep_sequences,1);
num_TMT = zeros(length_array,1);
for index = [1:length_array]
    num_TMT(index) = size(regexpi(pep_sequences{index},'K'),2)+1;
end
end

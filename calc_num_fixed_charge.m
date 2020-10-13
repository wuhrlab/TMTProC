function num_fixed_charge = calc_num_fixed_charge(pep_sequences)
length_array = size(pep_sequences,1);
num_fixed_charge = zeros(length_array,1);
for index = 1:length_array
     num_fixed_charge(index) = size(regexpi(pep_sequences{index},'[RKH]'),2)+1;
end
end
%function normalizes matrix of isotope envelope
%assumes that each envelope is in a row
function normalized_data=normalize_matrix_by_row(data)
normalized_data = data ./ repmat(sum(data, 2), 1, size(data,2)); 
end
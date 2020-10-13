% % this fucntion removes leading and trailing aa from peptides sequences and
% removes * (oxidized Methionin), replaces # with gg for GG-tag 
function cleaned_pep_sequences = clean_peptide_sequence(peptide_sequence_w_ends)

for index = 1:size(peptide_sequence_w_ends, 1)
    temp_seq = peptide_sequence_w_ends{index}; 
    temp_seq = temp_seq(3:(end-2));                 % get rid of leading "K." and  %get rid of trailing ".X" 
    temp_seq = strrep(temp_seq, '*','');            % get rid of * which codes for oxidized methionine
    temp_seq = strrep(temp_seq, ']','');            % get rid of N-terminal modification
    temp_seq = strrep(temp_seq,'#','');            % get rid of modification #
    temp_seq = strrep(temp_seq, 'X','i');           % Replace X (codes for I or L) with i 
    temp_seq = strrep(temp_seq, 'Z','e');           % Replace Z (codes for Q or E) with e
    temp_seq = strrep(temp_seq, 'B','d');           % Replace B (codes for D or N) with d
    cleaned_pep_sequences{index} = strrep(temp_seq, '#','gg');   % get rid of  # and replace with small gg
end
cleaned_pep_sequences = cleaned_pep_sequences';
end
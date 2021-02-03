function parameters=read_parameters(parameters_txt)


%read in file
str = fileread(parameters_txt);
%extract fields and values
pieces = regexp(str, '(?<field>\w+)=(?<value>[^\n]+)', 'names');
parameters = struct();  %variables

for k=1:length(pieces)
     field = pieces(k).field;
     val = pieces(k).value;
     
     tf= ismember(val, char([34 39]));    %remove " and ' from string
     val = strtrim(val(~tf));    %also remove spaces
     
     %try to convert to numeric.  If successful, use converted value,
     %otherwise use original (string) value
     valn = str2num(val);
     if ~isnan(valn) %conversion successful, use as numeric
        val = valn;
     end
     
     %add to output structure
     parameters.(field) = val;
end

end
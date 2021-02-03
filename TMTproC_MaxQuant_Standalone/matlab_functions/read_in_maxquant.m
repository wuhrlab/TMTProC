%  This function reads in search data from the MaxQuant output

function search_data = read_in_maxquant(filename) 

split_filename = regexp(filename,'\.','split');
if strcmp(split_filename{2},'txt')
    %get header and remove disallowed characters
    fileID = fopen(filename);
    header = fgetl(fileID);
    header = strsplit(strrep(strrep(strrep(strrep(header,' ','_'),')',''),'(',''),'/','_'),'\t');
    fclose(fileID);

    %read data
    search_data = tdfread(filename,'\t');

    %re-name headers
    search_data = cell2struct( struct2cell(search_data), header);
else
    fprintf('Data must be in .txt format')
    return
end

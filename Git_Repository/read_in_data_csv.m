%  This module reads in the data from Gygi's "CORE" format exported to csv
%  or xlsx/xls file
%  
%

function data = read_in_data_csv(filename) 


%read in excel file or csv file, xlsread is faster than reading in csv

split_filename = regexp(filename,'\.','split');
if strcmp(split_filename{2},'xls')|| strcmp(split_filename{2},'xlsx')
    [~,~,imported_data_cell_array] =  xlsread(filename);
elseif strcmp(split_filename{2},'csv')
    imported_data_cell_array = csvimport(filename);
    fprintf('Loaded Data from csv, excel format is faster to load \n')
else
    fprintf('Unknown file format. Can not load data. \n');
    return
end

save Imported_Data_Cell_Array
 
row_indeces= false(size(imported_data_cell_array,1),1);
row_indeces(2:end)=true; % indeces of imported_data.data has to be adjusted for missing header

col_indeces =0;
[~,col_indeces]             = ismember('ScanF',imported_data_cell_array(1,:));
if col_indeces
    data.ScanF                  = cell2mat(imported_data_cell_array(row_indeces, col_indeces));       % Read in ScanF => convert to number
else
    sprintf('Warning ScanF could not be found \n')
end


col_indeces =0;
[~,col_indeces]             = ismember('SrchID',imported_data_cell_array(1,:));
if col_indeces
    data.SearchID                  = cell2mat(imported_data_cell_array(2, col_indeces));       % Read in ScanF => convert to number
else
    sprintf('Warning SearchID could not be found \n')
end

col_indeces =0;
[~,col_indeces]             = ismember('Time',imported_data_cell_array(1,:));
if col_indeces
    data.Elution_Time                  = cell2mat(imported_data_cell_array(row_indeces, col_indeces));       % Read in ScanF => convert to number
else
    sprintf('Warning Elution Time could not be found \n')
end

col_indeces =0;
[~,col_indeces]             = ismember({'Precursor Intensity','PrecursorIntensity'}, imported_data_cell_array(1,:));
% get rid fo 0 entries
if sum(col_indeces)
    col_indeces(col_indeces==0)=[]; % remove 0 entries
    data.MS1_precursor_intensity            = cell2mat(imported_data_cell_array(row_indeces, col_indeces));       % Read in ScanF => convert to number
else
    sprintf('Warning MS1 Precursor intensity could not be found \n')
end

col_indeces =0;
[~,col_indeces]             = ismember({'Ion Injection Time','IonInjectionTime'},imported_data_cell_array(1,:));
if sum(col_indeces)
    col_indeces(col_indeces==0)=[];
    data.ion_injection_time            = cell2mat(imported_data_cell_array(row_indeces, col_indeces));       % Read in ScanF => convert to number
else
    sprintf('Warning Ion Injection Time could not be found \n')
end

%Calculate estimated ion numbers of precursors used for MS2 spectrum by
%denormalizing ion flux with ion injection time.
%data.denormalized_precursor_intensity = data.MS1_precursor_intensity.*data.ion_injection_time./1000;

col_indeces =0;
[~,col_indeces]             = ismember('Peptide',imported_data_cell_array(1,:));
if col_indeces
    data.pep_sequences          = imported_data_cell_array(row_indeces, col_indeces);          % Peptide sequence is red in
else
    sprintf('Warning Peptide Sequence could not be found \n')
end
col_indeces =0;

col_indeces =0;
[~,col_indeces]             = ismember('Reference',imported_data_cell_array(1,:));
if col_indeces
    data.prot_name          = imported_data_cell_array(row_indeces, col_indeces);          % Protein reference is red in
else
    sprintf('Warning Protein Reference could not be found \n')
end
col_indeces =0;


[~,col_indeces]             = ismember({'126 Sn','127 Sn','128 Sn','129 Sn','130 Sn','131 Sn','x179Sn', 'x180Sn'},imported_data_cell_array(1,:));
if sum(col_indeces)
    col_indeces(col_indeces ==0)=[];
    data.reporter_ions          = cell2mat(imported_data_cell_array(row_indeces, col_indeces));          % Read in S/N of rerporter Ions from data_data => convert to number 
else
    sprintf('Warning S/N levels for reporter ions could not be found \n')
end
col_indeces =0;
    
[~,col_indeces]             = ismember({'126 Mz Error','127 Mz Error','128 Mz Error','129 Mz Error','130 Mz Error','131 Mz Error'},imported_data_cell_array(1,:));
if col_indeces
    data.mz_err_reporter_ions   = cell2mat(imported_data_cell_array(row_indeces, col_indeces));       % mz_error of reporter Ions 
else
    sprintf('Warning m/z error for reporter ions could not be found \n')
end
col_indeces =0;

[~,col_indeces]             = ismember({'MPrecmin1 Sn','NPrec0 Sn','OPrecplus1 Sn','PPrecplus2 Sn','QPrecplus3 Sn','RPrecplus4 Sn','SPrecplus5 Sn','TPrecplus6 Sn','MPrecmin1Sn','NPrec0Sn','OPrecplus1Sn','PPrecplus2Sn','QPrecplus3Sn','RPrecplus4Sn','SPrecplus5Sn','TPrecplus6Sn';},imported_data_cell_array(1,:));
if sum(col_indeces)
    col_indeces(col_indeces==0)=[];
    data.precursor_iso_envelope = cell2mat(imported_data_cell_array(row_indeces, col_indeces));    % mz error of precursor Ions are 39:2:53
else
    sprintf('Precursor isotopic envelope could not be found \n')
end
col_indeces =0;
    
    
[~,col_indeces]             = ismember({'MPrecmin1 Mz Error','NPrec0 Mz Error','OPrecplus1 Mz Error','PPrecplus2 Mz Error','QPrecplus3 Mz Error','RPrecplus4 Mz Error','SPrecplus5 Mz Error','TPrecplus6 Mz Error';},imported_data_cell_array(1,:));
if col_indeces
    data.mz_err_precursor       = cell2mat(imported_data_cell_array(row_indeces, col_indeces));
else
    sprintf('Warning: Precursor isotopic envelope m/z error could not be found \n')
end
col_indeces =0;

[~,col_indeces]             = ismember({'AYnmin1Pmin1 Sn','BYnmin1P0 Sn','CYnmin1Pplus1 Sn','DYnmin1Pplus2 Sn','EYnmin1Pplus3 Sn','FYnmin1Pplus4 Sn','GYnmin1Pplus5 Sn','HYnmin1Pplus6 Sn','IYnmin1Pplus7 Sn','JYnmin1Pplus8 Sn','KYnmin1Pplus9 Sn','LYnmin1Pplus10 Sn','AYnmin1Pmin1Sn','BYnmin1P0Sn','CYnmin1Pplus1Sn','DYnmin1Pplus2Sn','EYnmin1Pplus3Sn','FYnmin1Pplus4Sn','GYnmin1Pplus5Sn','HYnmin1Pplus6Sn','IYnmin1Pplus7Sn','JYnmin1Pplus8Sn','KYnmin1Pplus9Sn','LYnmin1Pplus10Sn';},imported_data_cell_array(1,:));
if sum(col_indeces)
    col_indeces(col_indeces == 0)=[];
    data.Ynmin1_envelope        = cell2mat(imported_data_cell_array(row_indeces, col_indeces));    % Ynmin1 (cTMT) envelope
else
    sprintf('Warning: TMTc isotopic envelope could not be found \n')
end
col_indeces =0;

[~,col_indeces]             = ismember({'AYnmin1Pmin1 Mz Error','BYnmin1P0 Mz Error','CYnmin1Pplus1 Mz Error','DYnmin1Pplus2 Mz Error','EYnmin1Pplus3 Mz Error','FYnmin1Pplus4 Mz Error','GYnmin1Pplus5 Mz Error','HYnmin1Pplus6 Mz Error','IYnmin1Pplus7 Mz Error','JYnmin1Pplus8 Mz Error','KYnmin1Pplus9 Mz Error','LYnmin1Pplus10 Mz Error';},imported_data_cell_array(1,:));
if col_indeces
    data.mz_err_Ynmin1_envelope = cell2mat(imported_data_cell_array(row_indeces, col_indeces)); 
else
    sprintf('Warning: TMTc isotopic envelope m/z error could not be found \n')
end
col_indeces =0;

[~,col_indeces]             = ismember({'z'},imported_data_cell_array(1,:));
if col_indeces
    data.z                      = cell2mat(imported_data_cell_array(row_indeces, col_indeces));       % Charge
else
    sprintf('Warning: peptide charge state could not be found \n')
end
col_indeces =0;

[~,col_indeces]             = ismember({'PPM'},imported_data_cell_array(1,:));
if col_indeces
    data.ppm                    = cell2mat(imported_data_cell_array(row_indeces, col_indeces));  % ppm 
else
    sprintf('Warning: ppm could not be found \n')
end
col_indeces =0;

[~,col_indeces]             = ismember({'Theo m/z','TheoM_z'},imported_data_cell_array(1,:));
if sum(col_indeces)
    col_indeces(col_indeces ==0)=[];
    data.mz                     = cell2mat(imported_data_cell_array(row_indeces, col_indeces));      % m/z 
else
    sprintf('Warning: Theoretical m/z could not be found \n')
end
col_indeces =0;

[~,col_indeces]             = ismember({'Theo M+H','TheoM_H';},imported_data_cell_array(1,:));
if sum(col_indeces)
    col_indeces(col_indeces==0)=[];
    data.m_plus_H               = cell2mat(imported_data_cell_array(row_indeces, col_indeces));  % m+H 
else
    sprintf('Warning: Theoretical mass of peptide could not be found \n')
end
col_indeces =0;

[~,col_indeces]             = ismember({'PepID'},imported_data_cell_array(1,:));
if col_indeces
    data.PepID                  = cell2mat(imported_data_cell_array(row_indeces, col_indeces)); % Peptide ID
else
    sprintf('Warning: Peptide ID could not be found \n')
end
col_indeces =0;

[~,col_indeces]             = ismember({'SrchID'},imported_data_cell_array(1,:));
if col_indeces
    data.SrchID                 = cell2mat(imported_data_cell_array(row_indeces, col_indeces));   %Search ID
else
    sprintf('Warning: Search ID could not be found \n')
end
col_indeces =0;

% Convert m/z error for TMTc into approximate ppm error
data.TMTc_mz_approx = (data.m_plus_H-158.125818)./(data.z-1);

if exist('data.mz_err_Ynmin1_envelope','var')
    data.ppm_err_TMTc = data.mz_err_Ynmin1_envelope./(repmat(data.TMTc_mz_approx,1,size(data.mz_err_Ynmin1_envelope,2)));
end 
% From Search ID and Peptide ID generate a unique identifier for each Scan
% event
data.Unique_identifier = cell(length(data.PepID),1);
for index = 1: length(data.PepID)
   data.Unique_identifier{index} = strcat(num2str(data.SrchID(index)),'_',num2str(data.PepID(index)));
end

[~,col_indeces]             = ismember({'LDA Probability'},imported_data_cell_array(1,:));  %LDA Probability
if col_indeces
    data.LDA_probability        = cell2mat(imported_data_cell_array(row_indeces, col_indeces)); 
else
    sprintf('Warning: LDA probability could not be found \n')
end

col_indeces =0;

[~,col_indeces]             = ismember({'Isolation Mz','IsolationMz'},imported_data_cell_array(1,:));
if sum(col_indeces)
    col_indeces(col_indeces==0)=[];
    data.iso_mz             = cell2mat(imported_data_cell_array(row_indeces, col_indeces));  % mz value over which isolation window was centered
else
    sprintf('Warning: mz value for isolation window could not be found \n')
end

col_indeces =0;

[~,col_indeces]             = ismember({'fragment_sequence'},imported_data_cell_array(1,:));
if sum(col_indeces)
    col_indeces(col_indeces==0)=[];
    data.fragment_sequence = imported_data_cell_array(row_indeces, col_indeces);  % mz value over which isolation window was centered
else
    sprintf('Warning: fragment sequence could not be found \n')
end

[~,col_indeces]             = ismember({'missing_piece'},imported_data_cell_array(1,:));
if sum(col_indeces)
    col_indeces(col_indeces==0)=[];
    data.missing_piece = imported_data_cell_array(row_indeces, col_indeces);  % mz value over which isolation window was centered
else
    sprintf('Warning: missing piece could not be found \n')
end

[~,col_indeces]             = ismember({'num_TMT_fragment'},imported_data_cell_array(1,:));
if sum(col_indeces)
    col_indeces(col_indeces==0)=[];
    data.num_TMT_fragment = cell2mat(imported_data_cell_array(row_indeces, col_indeces));  % mz value over which isolation window was centered
else
    sprintf('Warning: num TMT fragment could not be found \n')
end

[~,col_indeces]             = ismember({'num_TMT_missing_piece'},imported_data_cell_array(1,:));
if sum(col_indeces)
    col_indeces(col_indeces==0)=[];
    data.num_TMT_missing_piece = cell2mat(imported_data_cell_array(row_indeces, col_indeces));  % mz value over which isolation window was centered
else
    sprintf('Warning: num TMT missing piece could not be found \n')
end


sprintf ('Data got loaded \n')
end
%  This module reads in the data from Gygi's "CORE" format exported to csv
%  or xlsx/xls file
%  
%

function data = read_in_data_csv(filename) 

%read in excel file or csv file, xlsread is faster than reading in csv

split_filename = regexp(filename,'\.','split');

fprintf('Loading data... \n')

if strcmp(split_filename{2},'xls')|| strcmp(split_filename{2},'xlsx')
    [~,~,imported_data_cell_array] =  xlsread(filename);
elseif strcmp(split_filename{2},'csv')
    imported_data_cell_array = csvimport(filename);
    fprintf('Loaded Data from csv, excel format is faster to load \n')
else
    fprintf('Unknown file format. Can not load data. \n');
    return
end

% save Imported_Data_Cell_Array

% Strip out any " from cell array in header line to address some
% compatibility issues with export from allosaurus

imported_data_cell_array (1,:) = cellfun(@(x) erase(x,'"'),imported_data_cell_array (1,:),'UniformOutput',false);
 
row_indeces= false(size(imported_data_cell_array,1),1);
row_indeces(2:end)=true; % indeces of imported_data.data has to be adjusted for missing header

col_indeces =0;
[~,col_indeces]             = ismember('MS/MS scan number',imported_data_cell_array(1,:));
if col_indeces
    data.ScanF                  = cell2mat(imported_data_cell_array(row_indeces, col_indeces));       % Read in ScanF => convert to number
else
    sprintf('Warning MS2 scan number could not be found \n')
end

col_indeces =0;
[~,col_indeces]             = ismember('Raw file',imported_data_cell_array(1,:));
if col_indeces
    data.SearchID                  = imported_data_cell_array(row_indeces, col_indeces);       % Read in ScanF => convert to number
else
    sprintf('Warning Raw file ID could not be found \n')
end

col_indeces =0;
[~,col_indeces]             = ismember('Retention time',imported_data_cell_array(1,:));
if col_indeces
    data.Elution_Time                  = cell2mat(imported_data_cell_array(row_indeces, col_indeces));       % Read in ScanF => convert to number
else
    sprintf('Warning Retention Time could not be found \n')
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


col_indeces =0;
[~,col_indeces]             = ismember('Sequence',imported_data_cell_array(1,:));
if col_indeces
    data.pep_sequences          = imported_data_cell_array(row_indeces, col_indeces);          % Peptide sequence is red in
else
    sprintf('Warning Peptide Sequence could not be found \n')
end
col_indeces =0;

col_indeces =0;
[~,col_indeces]             = ismember('Protein names',imported_data_cell_array(1,:));
if col_indeces
    data.prot_name          = imported_data_cell_array(row_indeces, col_indeces);          % Protein reference is red in
else
    sprintf('Warning Protein Names could not be found \n')
end
col_indeces =0;

%[~,col_indeces]             = ismember({'126 Sn','127n Sn','127c Sn','128n Sn','128c Sn','129n Sn','130n Sn', '130c Sn', '131n Sn','131c Sn','132n Sn','132c Sn','133n Sn','133c Sn',' 134n Sn'},imported_data_cell_array(1,:));
%if sum(col_indeces)
%    col_indeces(col_indeces ==0)=[];
%    data.reporter_ions          = cell2mat(imported_data_cell_array(row_indeces, col_indeces));          % Read in S/N of rerporter Ions from data_data => convert to number 
%else
%    sprintf('Warning S/N levels for reporter ions could not be found \n')
%end
%col_indeces =0;
    
%[~,col_indeces]             = ismember({'126 Mz Error','127n Mz Error','127c Mz Error','128n Mz Error','128c Mz Error','129n Mz Error','130n Mz Error', '130c Mz Error', '131n Mz Error','131c Mz Error','132n Mz Error','132c Mz Error','133n Mz Error','133c Mz Error',' 134n Mz Error'},imported_data_cell_array(1,:));
%if sum(col_indeces)
%    col_indeces(col_indeces ==0)=[];
%    data.mz_err_reporter_ions   = cell2mat(imported_data_cell_array(row_indeces, col_indeces));       % mz_error of reporter Ions 
%else
%    sprintf('Warning m/z error for reporter ions could not be found \n')
%end
%col_indeces =0;

[~,col_indeces]             = ismember({'mPrecmin1 Sn','nPrec0 Sn','oPrecplus1 Sn','pPrecplus2 Sn','qPrecplus3 Sn','rPrecplus4 Sn','sPrecplus5 Sn','tPrecplus6 Sn';},imported_data_cell_array(1,:));
if sum(col_indeces)
    col_indeces(col_indeces==0)=[];
    data.precursor_iso_envelope = cell2mat(imported_data_cell_array(row_indeces, col_indeces));    % mz error of precursor Ions are 39:2:53
else
    sprintf('Precursor isotopic envelope could not be found \n')
end
col_indeces =0;
    
    
[~,col_indeces]             = ismember({'mPrecmin1 Mz Error','nPrec0 Mz Error','oPrecplus1 Mz Error','pPrecplus2 Mz Error','qPrecplus3 Mz Error','rPrecplus4 Mz Error','sPrecplus5 Mz Error','tPrecplus6 Mz Error';},imported_data_cell_array(1,:));
if col_indeces
    data.mz_err_precursor       = cell2mat(imported_data_cell_array(row_indeces, col_indeces));
else
    sprintf('Warning: Precursor isotopic envelope m/z error could not be found \n')
end
col_indeces =0;

[~,col_indeces]             = ismember({'aYnmin1Pmin1 Sn','bYnmin1P0 Sn','cYnmin1Pplus1 Sn','dYnmin1Pplus2 Sn','eYnmin1Pplus3 Sn','fYnmin1Pplus4 Sn','gYnmin1Pplus5 Sn','hYnmin1Pplus6 Sn','iYnmin1Pplus7 Sn','jYnmin1Pplus8 Sn','kYnmin1Pplus9 Sn','lYnmin1Pplus10 Sn';},imported_data_cell_array(1,:));
if sum(col_indeces)
    col_indeces(col_indeces == 0)=[];
    data.Ynmin1_envelope        = cell2mat(imported_data_cell_array(row_indeces, col_indeces));    % Ynmin1 (cTMT) envelope
else
    sprintf('Warning: TMTc isotopic envelope could not be found \n')
end
col_indeces =0;

[~,col_indeces]             = ismember({'aYnmin1Pmin1 Mz Error','bYnmin1P0 Mz Error','cYnmin1Pplus1 Mz Error','dYnmin1Pplus2 Mz Error','eYnmin1Pplus3 Mz Error','fYnmin1Pplus4 Mz Error','gYnmin1Pplus5 Mz Error','hYnmin1Pplus6 Mz Error','iYnmin1Pplus7 Mz Error','jYnmin1Pplus8 Mz Error','kYnmin1Pplus9 Mz Error','lYnmin1Pplus10 Mz Error';},imported_data_cell_array(1,:));
if col_indeces
    data.mz_err_Ynmin1_envelope = cell2mat(imported_data_cell_array(row_indeces, col_indeces)); 
else
    sprintf('Warning: TMTc isotopic envelope m/z error could not be found \n')
end
col_indeces =0;

[~,col_indeces]             = ismember({'Charge'},imported_data_cell_array(1,:));
if col_indeces
    data.z                      = cell2mat(imported_data_cell_array(row_indeces, col_indeces));        % Charge
else
    sprintf('Warning: peptide charge state could not be found \n')
end
col_indeces =0;

%[~,col_indeces]             = ismember({'PPM'},imported_data_cell_array(1,:));
%if col_indeces
%    data.ppm                    = cell2mat(imported_data_cell_array(row_indeces, col_indeces));  % ppm 
%else
%    sprintf('Warning: ppm could not be found \n')
%end
%col_indeces =0;

[~,col_indeces]             = ismember({'Mass'},imported_data_cell_array(1,:));
if sum(col_indeces)
    col_indeces(col_indeces ==0)=[];
    data.m                     = cell2mat(imported_data_cell_array(row_indeces, col_indeces));      % m
    data.mz = data.m ./ data.z + 1.007276466812 ;
else
    sprintf('Warning: Theoretical m/z could not be found \n')
end
col_indeces =0;

%[~,col_indeces]             = ismember({'Theo M+H','TheoM_H';},imported_data_cell_array(1,:));
%if sum(col_indeces)
%    col_indeces(col_indeces==0)=[];
%    data.m_plus_H               = cell2mat(imported_data_cell_array(row_indeces, col_indeces));  % m+H 
%else
%    sprintf('Warning: Theoretical mass of peptide could not be found \n')
%end
%col_indeces =0;

[~,col_indeces]             = ismember({'Peptide ID'},imported_data_cell_array(1,:));
if col_indeces
    data.PepID                  = cell2mat(imported_data_cell_array(row_indeces, col_indeces)); % Peptide ID
else
    sprintf('Warning: Peptide ID could not be found \n')
end
col_indeces =0;

%[~,col_indeces]             = ismember({'SrchID'},imported_data_cell_array(1,:));
%if col_indeces
%    data.SrchID                 = cell2mat(imported_data_cell_array(row_indeces, col_indeces));   %Search ID
%else
%    sprintf('Warning: Search ID could not be found \n')
%end
%col_indeces =0;

[~,col_indeces]             = ismember({'Delta score'},imported_data_cell_array(1,:));  %Delta Score
if col_indeces
    data.delta_score        = cell2mat(imported_data_cell_array(row_indeces, col_indeces)); 
else
    sprintf('Warning: Delta score could not be found \n')
end

col_indeces =0;

[~,col_indeces]             = ismember({'Isolation Mz','IsolationMz'},imported_data_cell_array(1,:));
if sum(col_indeces)
    col_indeces(col_indeces==0)=[];
    data.iso_mz             = cell2mat(imported_data_cell_array(row_indeces, col_indeces));  % mz value over which isolation window was centered
else
    sprintf('Warning: mz value for isolation window could not be found \n')
end

% Convert m/z error for TMTc into approximate ppm error
data.TMTc_mz_approx = (data.m-158.125818)./(data.z-1);

data.ppm_err_TMTc = 1E6*data.mz_err_Ynmin1_envelope./(repmat(data.TMTc_mz_approx,1,size(data.mz_err_Ynmin1_envelope,2)));



% From Search ID and Peptide ID generate a unique identifier for each Scan
% event
data.Unique_identifier = cell(length(data.PepID),1);
for index = 1: length(data.PepID)
   data.Unique_identifier{index} = strcat(data.SearchID(index),'_',num2str(data.PepID(index)));
end


sprintf ('Data got loaded \n')
end
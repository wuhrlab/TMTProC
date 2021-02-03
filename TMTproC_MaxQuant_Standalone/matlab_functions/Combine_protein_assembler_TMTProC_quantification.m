%script links quantitative data from TMTProC with protein
%assembler
%TMTProC information is stored in exp structure, Protein assembler information
%is red in from excel or csv format


% read in data from protein assembler and TMTc structure
read_in_data = 1

Protein_assembler_filename = 'Fractionated_peptide_protein_matching.xlsx';


if read_in_data
    %load('combined_fractions_pep_export_2.mat'); %TMTc structure
    TMTProC_structure = results;
    split_filename = regexp(Protein_assembler_filename,'\.','split');
    if strcmp(split_filename{2},'xls')|| strcmp(split_filename{2},'xlsx')
        [~,~,pep_prot_map_cell_array] =  xlsread(Protein_assembler_filename);
    elseif strcmp(split_filename{2},'csv')
        pep_prot_map_cell_array = csvimport(Protein_assembler_filename);
        fprintf('Loaded Data from csv, excel format is faster to load \n')
    elseif strcmp(split_filename{2},'tsv')
        pep_prot_map_cell_array = dlmread(Protein_assembler_filename, '\t');
        fprintf('Loaded Data from tsv, excel format is faster to load \n')
    else
        fprintf('Unknown file format. Can not load data. \n');
        return
    end
end

%Find column in cell array from protein assembler which encodes peptide parsimony
[~,~,Pars_col_index]           = intersect({'Peptide_parsimony'},pep_prot_map_cell_array(1,:));
if isempty(Pars_col_index)
    sprintf('Warning Peptide Parsimony could not be found \n')
end

   %replace all NaN in peptide Parsimony with 'NaN' so that cell array function like ismember
   %can be exectuted
for index = 1: size(pep_prot_map_cell_array,1)
    if isnan(pep_prot_map_cell_array{index,Pars_col_index})
        pep_prot_map_cell_array{index,Pars_col_index} = 'NaN';
    end    
end



indeces_U_R_only = ismember(pep_prot_map_cell_array(:,Pars_col_index),{'R','U'});

%Clean out all peptides which are not U or R
%Remove all peptides which are not unique or razor in parsimony
pep_prot_map_cell_array_UR = pep_prot_map_cell_array(indeces_U_R_only,:);


%Find col_index Search ID in Protein Assembler
[~,~,SearchID_col_index]           = intersect({'SearchID'},pep_prot_map_cell_array(1,:));
if isempty(SearchID_col_index)
    sprintf('Warning SearchId could not be found \n')
end

%Find col index Peptide ID in Protein Assembler
[~,~,PeptideID_col_index]           = intersect({'PeptideID'},pep_prot_map_cell_array(1,:));
if isempty(PeptideID_col_index)
    sprintf('Warning PeptideId could not be found \n')
end

%Create unique peptide identifier by combining SearchID and Peptide ID for
%Protein Assembler
Unique_identifier = cell(length(pep_prot_map_cell_array_UR),1);
for index = 1: length(pep_prot_map_cell_array_UR)
   Unique_identifier{index} = strcat(num2str(pep_prot_map_cell_array_UR{index,SearchID_col_index}),'_',num2str(pep_prot_map_cell_array_UR{index,PeptideID_col_index}));
end


%Find colom for protein name in Protein Assembler
[~,~,Protein_col_index]           = intersect({'ProteinID'},pep_prot_map_cell_array(1,:));
if isempty(Protein_col_index)
    sprintf('Waring ProteinId could not be found \n')
end

%generate cell array with unique names of proteins
Unique_U_R_proteins = unique(pep_prot_map_cell_array_UR(:,Protein_col_index));

%Print out number of identified proteins
fprintf(strcat(sprintf('%d',size(Unique_U_R_proteins,1)),' Proteins were identified \n'))




% go through all unique proteins and calculate ratios
Quantified_Proteins = zeros(size(Unique_U_R_proteins,1),5);
Quantified_proteins_median = zeros(size(Unique_U_R_proteins,1),5);
Quantified_proteins_stdev = zeros(size(Unique_U_R_proteins,1),5);
How_many_peps_quantified = zeros(size(Unique_U_R_proteins,1),1);
%Only go through peptides in TMTc_structure which were actually quantified to speed up operation
quant_SN_TMTc = TMTProC_structure.ratios_times_SN(TMTProC_structure.indeces,:);
quant_unique_identifier_TMTc = TMTProC_structure.data.Unique_identifier(TMTProC_structure.indeces);

%Make Map Object with quant_unique_identifer_TMTc as key and
%quant_ratios_TMTc as content
mapObj_UniqueID_Quants = containers.Map(quant_unique_identifier_TMTc,num2cell(quant_SN_TMTc,2));


%Generate Map Object with Protein Name as key and the quant values as entry
mapObj_ProteinID_Quants = containers.Map;

for index = 1:length(pep_prot_map_cell_array_UR)
    %Check if MS spectra for unique Identifier was quantified => is in mapObj_UniqueID_Quants
    if isKey(mapObj_UniqueID_Quants,Unique_identifier{index}) 
        %Check if Protein ID already exists as key
       if isKey(mapObj_ProteinID_Quants,pep_prot_map_cell_array_UR{index,Protein_col_index})
           mapObj_ProteinID_Quants(pep_prot_map_cell_array_UR{index,Protein_col_index}) = [mapObj_ProteinID_Quants(pep_prot_map_cell_array_UR{index,Protein_col_index});mapObj_UniqueID_Quants(Unique_identifier{index})];
       %Create Protein ID as key with quantification as entry
       else
          mapObj_ProteinID_Quants(pep_prot_map_cell_array_UR{index,Protein_col_index}) = mapObj_UniqueID_Quants(Unique_identifier{index});
       end
    end
end

%Go through mapObj_ProteinID_Quants and calculate sum of SN for each protein and number of peptides 
for index = 1 : length(Unique_U_R_proteins)
    %Check if current protein is key to Unique_U_R_proteins (has been
    %quantified)
    if isKey(mapObj_ProteinID_Quants,Unique_U_R_proteins{index})
       %Calculate Mean
       Quantified_Proteins(index,1:10) = sum(mapObj_ProteinID_Quants(Unique_U_R_proteins{index}),1);
       How_many_peps_quantified(index) = size(mapObj_ProteinID_Quants(Unique_U_R_proteins{index}),1);
    end
    index
end
output = [Unique_U_R_proteins,num2cell(Quantified_Proteins),num2cell(How_many_peps_quantified)];
output = [{'Protein ID','TMTPro0','126','128c','129n','130c','131n','131c','133c','134n','HeavyTMTPro','# quant peps'};output];
xlswrite('Protein_ID_Quantified_Proteins_4',output)

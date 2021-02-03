function TMTproC_maxquant_compatible(input_file)

%addpath(genpath('matlab_functions'))

fprintf('Beginning deconvolution \n')
params = read_parameters('parameters/tmtproc_input_params.txt');

% Load Isolation window
load('Saved_Iso_Window_0p4')

%load TMT isotopic Impurity Matrix
load('TMTPro_Impurity_Matrix_2019')

%Use_Precursor = 0;
%which_channels_to_use = [0,1,1,1,1,1,1,1,1,0]; % binary which channels to use, first position is TMTPro0, second position is TMTPro126, 9th position is TMTPro134N, 10th position is currently not usable.
%Noiseband = 1; % Conversion from S:N to ions. Dependent on the resolution used
%SN_cutoff = 0; 
%Use_ppm_filter = 0;  %Use ppm filter as defined by Alex Johnson => throw out all spectra in which at least one compliment ion differs by more than 10ppm from expected values


%Call the master program
for filenames = string(input_file)
    results = demo_TMTProc_iso_Window(filenames{1},params.noiseband,params.which_channels_to_use,params.sn_cutoff,Array_Iso_Window',params.use_ppm_filter,TMT_Impurity_Matrix,params.use_precursor);
    saving_results_name = strrep(strcat(strsplit(filenames{1},'.'),'_results'),'_intermediate','')
    save(saving_results_name{1},'results')
    
    data_table = struct2table(results.data);
    
    to_export = horzcat(results.ratios(:,1:10),results.data.sum_ions_Ynmin1, results.passed_ppm_filter, results.indeces, results.mobile_protons_index);
    varnames ={'ratios0','ratio1','ratio2','ratio3','ratio4','ratio5','ratio6','ratio7','ratio8','ratios9','sum_SN','passed_ppm_filter','indeces','mobile_protin_index'};
    for i=1:length(varnames)
        data_table.(varnames{i}) = to_export(:,i);
    end
    
    writetable(data_table, strcat(saving_results_name{1},'.csv'))
    
end

end
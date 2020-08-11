% Define the isolation window for MS2 spectra

% Load Isolatino window
% Use load('Saved_Iso_Window_0p5') for 0.5 Isolation window on Lumos (TGR)
% Use load('Saved_Iso_Window_0p6') for 0.6 Isolation window on Lumos (TGR)
% Use load('Saved_Iso_Window_1p0') for 1.0 Isolation window on Lumos (TGR)
load('Saved_Iso_Window_0p4')
%load TMT isotopic Impurity Matrix
% Use load('TMT_Impurity_Matrix_TMT0.mat') for TMT0
% Use load('TMT_Impurity_Matrix_2012.mat') for isotopic impurities
% collected in 2012
% Use load('TMT_Impurity_Matrix_2016.mat') for isotopic impurities
% collected in 2016
load('TMTPro_Impurity_Matrix_2019')
Use_Precursor = 0
which_channels_to_use = [0,1,1,1,1,1,1,1,1,0]; % binary which channels to use, first position is TMTPro0, second position is TMTPro126, 9th position is TMTPro134N, 10th position is currently not usable. 

%Call the master program
for filenames = {'TGR_09108_HeLa_YeastInterference_1to1_2uL_TMTc_noFAIMS.xlsx', 'TGR_09097_HeLa_YeastInterference_1to1_2uL_TMTc_FAIMS_-45V_-65V.xlsx'}
    results = demo_TMTProc_iso_Window(filenames{1},1,which_channels_to_use,1,100,Array_Iso_Window',0.1,TMT_Impurity_Matrix,Use_Precursor);
    saving_results_name = strsplit(filenames{1},'.')
    save(saving_results_name{1},'results')
    clear results
end

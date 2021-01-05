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
% Use load('TMT_Impurity_Matrix_2019.mat') for isotopic impurities
% collected in 2019

load('TMTPro_Impurity_Matrix_2019.mat')

Use_Precursor = 0;
which_channels_to_use = [0,1,1,1,1,1,1,1,1,0]; % binary which channels to use, first position is TMTPro0, second position is TMTPro126, 9th position is TMTPro134N, 10th position is currently not usable.
Noiseband = 1;
SN_cutoff = 0;
Use_ppm_filter = 0;  %Use ppm filter as defined by Alex Johnson => throw out all spectra in which at least one compliment ion differs by more than 10ppm from expected values

%Call the master program
for filenames = {'scan11131.xlsx'}
    results = demo_TMTProc_iso_Window(filenames{1},Noiseband,which_channels_to_use,SN_cutoff,Array_Iso_Window',Use_ppm_filter,TMT_Impurity_Matrix,Use_Precursor);
    saving_results_name = strsplit(filenames{1},'.')
    save(saving_results_name{1},'results')
    % plot summary of results
    plot_run_summaries(results);
end


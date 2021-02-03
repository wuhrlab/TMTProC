% results = demo_TMTProc_iso_Window(filenames{1},1,which_channels_to_use,Noiseband,SN_cutoff,Array_Iso_Window',0.1,TMT_Impurity_Matrix,Use_Precursor)
% Master Function to run TMTProC quantification 
% 
%  Martin  Wuehr,  2012-2020 
% 
function exp=demo_TMTProc_iso_Window(filename,noiseband,which_channels_used,sum_ions_Ynmin1_cutoff,Iso_Window,use_ppm_filter,TMT_Impurity_Matrix,Use_Precursor) 

% define filename
identifier = regexp(filename, '\.', 'split');       % convert filename to identifier 

exp.Iso_Window = Iso_Window;
exp.use_ppm_filter = use_ppm_filter;
exp.which_channels_used = which_channels_used; %Which channels are used logic 10 element vector [TMTPro0,126,127,...] 
exp.noiseband = noiseband; %How many charges is the noise band? for elite 5 (15k), 3.5 (30k), 2.5 (60k), for QExactive 5(7.5k), 3.5(15k), 2(30k)
exp.sum_ions_Ynmin1_cutoff = sum_ions_Ynmin1_cutoff; %Number of ions in TMTc cluster


exp.data =  read_in_data_csv(filename);         % read in data if not defined yet, 
exp.name =  identifier{1}; 

%exp.data.pep_sequences = clean_peptide_sequence(exp.data.pep_sequences);  % clean up sequeences (remove dots and trailing and leading AAs) 
exp.data.num_TMT = calculate_num_TMT(exp.data.pep_sequences);             % count number of Ks to derive number of TMT labels
exp.data.num_fix_charge = calc_num_fixed_charge(exp.data.pep_sequences);  % count number of N-term, K, R and H to calculate fixed charges

if exist('exp.data.precursor_iso_envelope','var')
    exp.data.sum_precursor_envelope_norm = sum(exp.data.precursor_iso_envelope, 2) ./ exp.data.z;   % Calculate sum of envelope for precursor
end
%exp.data.sum_Ynmin1_envelope_norm = sum(exp.data.Ynmin1_envelope, 2) ./ (exp.data.z - 1);        % Calculate sum of envelope for Y(n-1) 

%exp.data.reporter_ions_normalized = normalize_matrix_by_row(exp.data.reporter_ions);        % normalize Reporter Ions to sum of reporter ions
exp.num_seqs = size(exp.data.pep_sequences, 1);          % number of peptides 
% Calculate theoretical precursor envelope:
exp.theoretical_precursor_envelope = zeros(exp.num_seqs, 8);  
%show waitbar for progress of isotopic precursor cla
h=waitbar(0,'Calculating peptide isotopic envelope');
% load previously calculated isotopic envelopes
if exist('matlab_functions/saved_isotopic_envelopes.mat','file')
   load('matlab_functions/saved_isotopic_envelopes.mat')
   fprintf('Loaded saved isotopic envelopes \n')
else
   Previously_saved_IsoEnvs = containers.Map; 
end

for index = 1:exp.num_seqs
    %Check if Isotopic envelope for peptide sequence has already been
    %calculated
    if isKey(Previously_saved_IsoEnvs,exp.data.pep_sequences{index})
        exp.data.theoretical_precursor_envelope(index,:) = Previously_saved_IsoEnvs(exp.data.pep_sequences{index});
        % If isotopic envelope has not been calculated yet calculate now
    else
        MD = isotopicdist(exp.data.pep_sequences(index), 'fftresolution', 20, 'SHOWPLOT', false);    % Calculate isotope envelope of the peptide alone
        MD = [MD(:,2); zeros(8,1)];                 % Make sure that theoretical isotope vector contains at least 12 entries
        exp.data.theoretical_precursor_envelope(index,:) = normalize_matrix_by_row(MD(1:12)'); % normalize envelope to sum=1
        % Store the newly calculated isotopic envelope in the map object
        % for future use
        Previously_saved_IsoEnvs(exp.data.pep_sequences{index}) = exp.data.theoretical_precursor_envelope(index,:);
    end
    
    waitbar(index/exp.num_seqs);
end
close(h);

%save('matlab_functions/saved_isotopic_envelopes2.mat','Previously_saved_IsoEnvs'); 


[exp.ratios,exp.goodnes_fit,exp.data.sum_ions_Ynmin1,exp.calculated_TMTc_isotope_pattern,exp.Percentage_Permeates_at_Position,exp.Precursor_Matrices,exp.binary_abundand_Ynmin1_positions, exp.precursor_used_for_isolation_window] = main_function_w_iso_window_TMTPro(exp.data,exp.which_channels_used,exp.Iso_Window,exp.noiseband,TMT_Impurity_Matrix,Use_Precursor);

% multiply ratios by signal input
exp.ratios_times_SN = exp.ratios.*exp.data.sum_ions_Ynmin1;


% calculate theoretical precursor distribution (vector) for all calculated mixing ratios
%exp.theoretical_precursor_w_TMT = batch_calculate_theoretical_precursor_envelope_TMT(exp);

%Filter out peptides that do not pass ppm filter for comp ion
%quantificaiton
if exp.use_ppm_filter
    exp = does_ppm_filter_pass(exp);
else
    exp.passed_ppm_filter = true(length(exp.data.z),1);
end

% Calculate mobile protons
exp.data.mobile_protons = exp.data.z - exp.data.num_fix_charge;
exp.mobile_protons_index = exp.data.mobile_protons >= 1;

% Calculate which peptide pass quality criteria
exp.indeces = exp.passed_ppm_filter & exp.data.sum_ions_Ynmin1 > exp.sum_ions_Ynmin1_cutoff & exp.data.mobile_protons < 1;



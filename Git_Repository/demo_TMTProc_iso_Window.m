% demo_TMTc_iso_Window('v11104_Ynmin1_recall.csv',3.5,[1,1,1,1,1],1,1000,Array_Iso_Window,0.0005,1)
%demo_TMTc_iso_Window(filename,noiseband,which_channels,calculate_cosine_distance,minimum
%number ions in TMTc envelope,goodnesfit_cutoff, Function_Iso_Window,use  precursor)
% Master Function to run Ynmin quantification 
% demo_Ynmin1 takes in filename which contains quantitative data exported
% from CoreQuant and the noiseband. The noiseband is how many charges
% equivalent the noiseband for the experiment is. E.g. OrbitrapElite 30k
% (400mz) experiment => noiseband=3.5, 
%
%  Martin  Wuehr,  2012 
% 
function exp=demo_TMTProc_iso_Window(filename,noiseband,which_channels_used,calculate_cosine_distance,sum_ions_Ynmin1_cutoff,Iso_Window,goodnesfit_cutoff,TMT_Impurity_Matrix,Use_Precursor) 

% define filename
identifier = regexp(filename, '\.', 'split');       % convert filename to identifier 

exp.Iso_Window = Iso_Window;
exp.which_channels_used = which_channels_used; %Which channels are used logic 10 element vector [126,127,128,129/130,131] 
exp.noiseband = noiseband; %How many charges is the noise band? for elite 5 (15k), 3.5 (30k), 2.5 (60k), for QExactive 5(7.5k), 3.5(15k), 2(30k)
exp.sum_ions_Ynmin1_cutoff = sum_ions_Ynmin1_cutoff; %Number of ions in TMTc cluster
exp.calculate_cosine_distance = calculate_cosine_distance;              % Calculate cosine distance for samples with known mixing ratios
exp.known_ratios = [0,2.19966067289355,2.66815201102505,2.97710822343980,1.90166524696685,2.25207055290412,2.71685072821662,2.44227988560770,2.74160519047206,0];    % Known mixing ratios when analyzing standards
exp.data =  read_in_data_csv(filename);         % read in data if not defined yet, 
exp.name =  identifier{1}; 
exp.data.pep_sequences = clean_peptide_sequence(exp.data.pep_sequences);  % clean up sequeences (remove dots and trailing and leading AAs) 
exp.data.num_TMT = calculate_num_TMT(exp.data.pep_sequences);             % count number of Ks to derive number of TMT labels
exp.data.num_fix_charge = calc_num_fixed_charge(exp.data.pep_sequences);  % count number of N-term, K, R and H to calculate fixed charges
if exist('exp.data.precursor_iso_envelope','var')
    exp.data.sum_precursor_envelope_norm = sum(exp.data.precursor_iso_envelope, 2) ./ exp.data.z;   % Calculate sum of envelope for precursor
end
%exp.data.sum_Ynmin1_envelope_norm = sum(exp.data.Ynmin1_envelope, 2) ./ (exp.data.z - 1);        % Calculate sum of envelope for Y(n-1) 
exp.data.reporter_ions_normalized = normalize_matrix_by_row(exp.data.reporter_ions);        % normalize Reporter Ions to sum of reporter ions
exp.num_seqs = size(exp.data.pep_sequences, 1);          % number of peptides 
% Calculate theoretical precursor envelope:
exp.theoretical_precursor_envelope = zeros(exp.num_seqs, 8);  
%show waitbar for progress of isotopic precursor cla
h=waitbar(0,'Calculating peptide isotopic envelope');
for index = 1:exp.num_seqs 
    MD = isotopicdist(exp.data.pep_sequences(index), 'fftresolution', 20, 'SHOWPLOT', false);    % Calculate isotope envelope of the peptide alone 
    MD = [MD(:,2); zeros(8,1)];                 % Make sure that theoretical isotope vector contains at least 12 entries
    exp.data.theoretical_precursor_envelope(index,:) = normalize_matrix_by_row(MD(1:12)');   % normalize envelope to sum=1
    waitbar(index/exp.num_seqs);
end
close(h);

save temp; 

[exp.ratios,exp.goodnes_fit,exp.data.sum_ions_Ynmin1,exp.calculated_TMTc_isotope_pattern,exp.Percentage_Permeates_at_Position,exp.Precursor_Matrices,exp.binary_abundand_Ynmin1_positions, exp.precursor_used_for_isolation_window] = main_function_w_iso_window_TMTPro(exp.data,exp.which_channels_used,exp.Iso_Window,exp.noiseband,TMT_Impurity_Matrix,Use_Precursor);
% multiply ratios with 20
exp.ratios=exp.ratios.*20;
save temp;
if exp.calculate_cosine_distance                 % calculate cosine distance for know mixing ratios  
    exp.cosine_distance = pep_cosine_distance(exp.ratios,exp.known_ratios);
end

% calculate theoretical precursor distribution (vector) for all calculated mixing ratios
%exp.theoretical_precursor_w_TMT = batch_calculate_theoretical_precursor_envelope_TMT(exp);


if exist('goodnesfit_cutoff','var')
    exp.goodnesfit_cutoff = goodnesfit_cutoff; %Check if goodnesfit cutoff is supplied otherwise calculate cutoff based on minimum TMTc signal
else
    %Based on labbook entry 2012-10-05 or Suppl Fig. 4F in TMTc paper 
    %Diff < a*sum(Ions in TMTc)^b
    exp.a_quality_cutoff = 1.9188; %m for cutoff criterion
    exp.b_quality_cutoff = -0.8614; %t for cutoff criterion
    exp.goodnesfit_cutoff = exp.a_quality_cutoff.*(exp.sum_ions_Ynmin1_cutoff.^exp.b_quality_cutoff);
end

% Calculate mobile protons
exp.data.mobile_protons = exp.data.z - exp.data.num_fix_charge;
exp.mobile_protons_index = exp.data.mobile_protons >= 1;

% Calculate which peptide pass quality criteria
exp.indeces = exp.goodnes_fit(:,1) <exp.goodnesfit_cutoff & exp.data.sum_ions_Ynmin1 > exp.sum_ions_Ynmin1_cutoff;% & exp.data.mobile_protons < 1;

plot_only_yeast = 0;
if plot_only_yeast %Plot only ratios for yeast peptides
    restrictive_string = 'ORF';
    exp.indeces_yeast_only = strncmpi(exp.data.prot_name, restrictive_string, 3);
    exp.indeces = exp.indeces & exp.indeces_yeast_only;
end


%plot run summaries

exp = plot_run_summaries(exp); 

% any2csv(data,',',0,strcat('results_',identifier{1},'.csv'));       %  write structure to csv file


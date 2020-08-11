%
% this function goes through data structure containing precurser envelope (8) and
% Y(n-1)_envelope (12) and returns ratios (5) and goodness of fit(1)
%
function [ratios,goodnes_fit,sum_ions_Ynmin1,calculated_TMTc_isotope_pattern,Percentage_Permeates_at_Position_stored,Precursor_Matrices,binary_abundand_Ynmin1_positions_saved, precursor_used_for_isolation_window] = main_function_w_iso_window_TMTPro(data,which_channels_used,Array_Iso_Window,noiseband,TMT_impurity_Matrix,Use_Precursor)
%Normalize each Y(n-1) envelope to 1
Ynmin1_envelope_norm = [normalize_matrix_by_row(data.Ynmin1_envelope(:,:)),zeros(size(data.Ynmin1_envelope,1),3)];
% How many for loops to perform?
max_for = size(data.ScanF, 1);
%initialize ratios and goodness_fit
ratios=zeros(max_for,10);
goodnes_fit=zeros(max_for,1);
calculated_TMTc_isotope_pattern = zeros(max_for,15);

%initiate_Storage_for_binary_abundance
binary_abundand_Ynmin1_positions_saved = zeros(max_for,15);

%initiate What percentage of each position permeates the isolation window
Percentage_Permeates_at_Position_stored = zeros(max_for,12);

%Initialize structure to store the TMT postion matrix P(TMT) P_126...P_121
temp = cell(max_for,1);
Precursor_Matrices = struct('P_126',temp,'P_127',temp,'P_128',temp,'P_130',temp,'P_131',temp);

% initializing Sum Ynmin 1
sum_ions_Ynmin1 = zeros(length(data.ScanF),1);
% initializes waitbar
%h = waitbar(0,'Calculating ratios based on TMTc cluster');
% replaced waitbar with parfor progress
parfor_progress(max_for); % Initialize

% copy out data from structure so that it can be sliced for forloop
%Calculate what percentage of each isotope makes it through the
%isotpic envelope

if isfield(data,'iso_mz')
    offset = data.iso_mz-data.mz; % assume that the center of the isolation window is placed on the pseudo-monoisotopic peak
else
        offset = zeros(length(data.z));
end

z = data.z;
Ynmin1_envelope = data.Ynmin1_envelope;
measured_precursor_envelope = data.precursor_iso_envelope;
theoretical_precursor_envelope = data.theoretical_precursor_envelope;
num_TMT = data.num_TMT;

%Initialize vector to store if isolation window was used. Keep outside data
%structure as it has to be updated during parfor loop. 
precursor_used_for_isolation_window = zeros(length(data.num_TMT),1);

for index = 1:max_for      % Can be switched between parfor and for loop for speed/debugging
    if index == 0
        1; %Debugging possibility to stop
    end
    parfor_progress; % Count
    Position_Precursor = [-1:10]/z(index)-offset(index); %adjusts the isolation window for the different charge states
    Which_mz_values_from_Iso_Window = zeros(size(Position_Precursor));
    for index_2 = 1:length(Position_Precursor)
        [~,Which_mz_values_from_Iso_Window(index_2)] = min(abs(Position_Precursor(index_2)-Array_Iso_Window(:,1)));
    end
    Percentage_Permeates_at_Position = Array_Iso_Window(Which_mz_values_from_Iso_Window,2);
    Percentage_Permeates_at_Position = Percentage_Permeates_at_Position./max(Percentage_Permeates_at_Position);
    
    % Check if enough precursor was measured, if so and if precursor flag
    % has been selected modify the permeates position by what was measured
    if Use_Precursor && sum(measured_precursor_envelope(index,:)) > 50
       %Update Array if precursor was used to define isolation window
       precursor_used_for_isolation_window(index) = 1;
       Which_precursor_positions_to_consider = logical(Percentage_Permeates_at_Position);
       %Find highest and lowest position that can be populated based on
       %isolation window => add one more position to it
       old_first = find(Which_precursor_positions_to_consider, 1 ); 
       new_first = old_first -1; 
       if new_first < 1
           new_first = 1;
       end
       Which_precursor_positions_to_consider(new_first) = 1;
       
       old_last = find(Which_precursor_positions_to_consider,1,'last'); 
       new_last = old_last + 1; 
       if new_last > 12
           new_last = 12;
       end
       Which_precursor_positions_to_consider(new_last) = 1;
       % Pull out surving precursor
       Percentage_Permeates_at_Position = [measured_precursor_envelope(index,:),0,0,0,0]'.* Which_precursor_positions_to_consider;
       Percentage_Permeates_at_Position = Percentage_Permeates_at_Position./sum(Percentage_Permeates_at_Position);
       
    end
    
    Percentage_Permeates_at_Position_stored(index,:) = Percentage_Permeates_at_Position';
    
    %Calculate each TMT postion matrix P(TMT) for given isotopic
    %precursor envelope and number of TMTs labeled with considering
    %charge state and isolation window
    [P_0,P_1,P_2,P_3,P_4,P_5,P_6,P_7,P_8,P_9] = calc_all_TMTPro_position_and_isotope_distribution_w_iso_W(theoretical_precursor_envelope(index,:), num_TMT(index),Percentage_Permeates_at_Position,TMT_impurity_Matrix);
    %Store the calculated values
    Precursor_Matrices(index).P_0 = P_0;
    Precursor_Matrices(index).P_1 = P_1;
    Precursor_Matrices(index).P_2 = P_2;
    Precursor_Matrices(index).P_3 = P_3;
    Precursor_Matrices(index).P_4 = P_4;
    Precursor_Matrices(index).P_5 = P_5;
    Precursor_Matrices(index).P_6 = P_6;
    Precursor_Matrices(index).P_7 = P_7;
    Precursor_Matrices(index).P_8 = P_8;
    Precursor_Matrices(index).P_9 = P_9;
    
    % Calculate Ynmin1 vector for 0.2:0.2:0.2:0.2:0.2 mixing ratios. Use
    % only positions which have more than 1% abundance for fitting later
    % Use binary that covers entire 15 vector
    binary_vector_all=ones(1,15);
    ratios_equally_distributed = which_channels_used./sum(which_channels_used); %Calculate the vector for the equally distributed example for a given number of channels
    binary_abundand_Ynmin1_positions = 0.01 < calculate_theoretical_Ynmin1_from_matrix(ratios_equally_distributed,P_0,P_1,P_2,P_3,P_4,P_5,P_6,P_7,P_8,P_9,binary_vector_all);
    % Save for later
    binary_abundand_Ynmin1_positions_saved(index,:) = binary_abundand_Ynmin1_positions; 
    
    %update data.sum_ions_Ynmin1 considering only the positions which
    %would be filled by a 1:1:1 envelope
    sum_ions_Ynmin1(index) = noiseband*sum([Ynmin1_envelope(index,:),0,0,0].*binary_abundand_Ynmin1_positions);
    %exp.data.sum_ions_Ynmin1=sum(data.Ynmin1_envelope(:,2:8), 2) ./ (data.z - 1).*exp.; % The charges rather than ions version: exp.data.sum_ions_Ynmin1=sum(exp.data.Ynmin1_envelope, 2)./exp.noiseband;  %Number of charges in Ynmin1 envelope % only consider positions 0 to 10 as sometimes the 12th position will be cutof => 0 charges

    if sum_ions_Ynmin1(index) > 50
        %renormalize Ynmin1 envelope only for positions considered
        Ynmin1_envelope_renormalized=normalize_matrix_by_row(binary_abundand_Ynmin1_positions.*Ynmin1_envelope_norm(index,:));
        %options defines the tolerance for the solutions (default would be e-6)
        options = optimset('TolFun', 1e-7, 'TolX', 1e-6,'Display', 'off');
        %difference is the function to optimize (square distance between
        %observed and calculated Y(n-1)
        difference = @(x)sum(binary_abundand_Ynmin1_positions.*((Ynmin1_envelope_renormalized-calculate_theoretical_Ynmin1_from_matrix(x,P_0,P_1,P_2,P_3,P_4,P_5,P_6,P_7,P_8,P_9,binary_abundand_Ynmin1_positions)).^2));
        
        %problem definitions: 'x0' possible starting point, 'objective' defines the
        %function that will be minimized.'lb' lower bound for x,'ub' upper bound
        %for x, 'Aeq' and 'beq' are defining the requirement that the sum of all
        %five x [1 1 1 1 1] is 1, options handles the tolerance for the solver.
        problem = createOptimProblem('fmincon','x0', ratios_equally_distributed,'objective',difference,'lb',[0,0,0,0,0,0,0,0,0,0],'ub',which_channels_used,'Aeq',[1,1,1,1,1,1,1,1,1,1], 'beq',1, 'options',options); %defines problem
        [ratios(index,:), goodnes_fit(index)] = fmincon(problem);
        calculated_TMTc_isotope_pattern(index,:) = calculate_theoretical_Ynmin1_from_matrix(ratios(index,:),P_0,P_1,P_2,P_3,P_4,P_5,P_6,P_7,P_8,P_9,binary_abundand_Ynmin1_positions);
    end
end

parfor_progress(0); % Close parfor bar
end
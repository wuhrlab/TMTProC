%This function takes in the mixing ratio of 5plex TMT (z) (126,127,128,130,131), the
%precursor_isotopic_distribution_wo_TMT (length=8)
%for the peptide and the number of TMT (1) it is labeled with,
%num_TMT. It calculates a Matrix that encodes Isotopic envelope position
%starting from pseudo-Monoistopic -1 to pseudo-Mono +6
%The input of z should add up to 1 output will add up to 1

function [P_0,P_1,P_2,P_3,P_4,P_5,P_6,P_7,P_8,P_9]= calc_all_TMTPro_position_and_isotope_distribution_w_iso_W(isotope_prec_wo_TMT, num_TMT,Percentage_Permeates_at_Position,TMT_Impurity_Matrix)



Matrix_Percentage_Permeates_at_Position = repmat(Percentage_Permeates_at_Position',11,1);
%For each of the TMT-channels calculate the precursor TMT matrix P(126 ...
%131) and simulate the isolation window for the MS2 spectrum
P_0 = convolve_precursor_w_TMT(isotope_prec_wo_TMT,TMT_Impurity_Matrix.Impurities_L0,num_TMT);
P_0 = P_0 .* Matrix_Percentage_Permeates_at_Position;
P_1 = convolve_precursor_w_TMT(isotope_prec_wo_TMT,TMT_Impurity_Matrix.Impurities_L1,num_TMT);
P_1 = P_1 .* Matrix_Percentage_Permeates_at_Position;
P_2 = convolve_precursor_w_TMT(isotope_prec_wo_TMT,TMT_Impurity_Matrix.Impurities_L2,num_TMT);
P_2 = P_2 .* Matrix_Percentage_Permeates_at_Position;
P_3 = convolve_precursor_w_TMT(isotope_prec_wo_TMT,TMT_Impurity_Matrix.Impurities_L3,num_TMT);
P_3 = P_3 .* Matrix_Percentage_Permeates_at_Position;
P_4 = convolve_precursor_w_TMT(isotope_prec_wo_TMT,TMT_Impurity_Matrix.Impurities_L4,num_TMT);
P_4 = P_4 .* Matrix_Percentage_Permeates_at_Position;
P_5 = convolve_precursor_w_TMT(isotope_prec_wo_TMT,TMT_Impurity_Matrix.Impurities_L5,num_TMT);
P_5 = P_5 .* Matrix_Percentage_Permeates_at_Position;
P_6 = convolve_precursor_w_TMT(isotope_prec_wo_TMT,TMT_Impurity_Matrix.Impurities_L6,num_TMT);
P_6 = P_6 .* Matrix_Percentage_Permeates_at_Position;
P_7 = convolve_precursor_w_TMT(isotope_prec_wo_TMT,TMT_Impurity_Matrix.Impurities_L7,num_TMT);
P_7 = P_7 .* Matrix_Percentage_Permeates_at_Position;
P_8 = convolve_precursor_w_TMT(isotope_prec_wo_TMT,TMT_Impurity_Matrix.Impurities_L8,num_TMT);
P_8 = P_8 .* Matrix_Percentage_Permeates_at_Position;
P_9 = convolve_precursor_w_TMT(isotope_prec_wo_TMT,TMT_Impurity_Matrix.Impurities_L9,num_TMT);
P_9 = P_9 .* Matrix_Percentage_Permeates_at_Position;



%For a given mixing ratio x calculate the Mixing ratio Matrix P_M as a
%weighted sum of P_126 to P_131
%Precursor_w_TMT_Matrix_for_mixing_ratio=x(1).*P_126+x(2).*P_127+x(3).*P_128+x(4).*P_130+x(5).*P_131;
% if 0
%     bar(Matrix_precursor_TMT_isotopes')
%     legend('delta150','delta151','delta152','delta153','delta154', 'delta155','delta156');
%     xlabel('Position in Isotopic Envelope starting at -1 to +6')
% end


		
		
		

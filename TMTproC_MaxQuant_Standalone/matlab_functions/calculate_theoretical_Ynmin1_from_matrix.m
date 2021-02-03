% This function calculates the theoretical Y_nmin1 based on a given mixing
% ratio (x,5), the Precursor_matrix for each TMT channel P(126 to 131). P,
% which is calculated based on AA sequence and number of TMTs labeled with,
% has 6 rows coding for delta masses, 12 columns coding for the position in
% the precursor isotopic envelope (from -1 to +10 with the pseudo
% mono-isotopic peak defining position 0)
% Function returns Theoretical_Ynmin1 envelope (vector of 15 from 0
% (Monoisotopic peak (0) to +14) and a binary vector of same length
% encoding which positions 

function theoretical_Ynmin1=calculate_theoretical_Ynmin1_from_matrix(x,P_0,P_1,P_2,P_3,P_4,P_5,P_6,P_7,P_8,P_9,binary_abundand_Ynmin1_positions)

%For a given mixing ratio x calculate the Mixing ratio Matrix P_M as a
%weighted sum of P_126 to P_131
Precursor_w_TMT_Matrix_for_mixing_ratio=x(1).*P_0+x(2).*P_1+x(3).*P_2+x(4).*P_3+x(5).*P_4+x(6).*P_5+x(7).*P_6+x(8).*P_7+x(9).*P_8+x(10).*P_9;

% Go through the matrix and sum up vector Y(n-1), position coded by k,
% which fullfils the condition i+k+2=j
theoretical_Ynmin1=zeros(1,15);
for k = 1:15
    for i=1:11
        for j=1:12
            if i+k-9==j
                theoretical_Ynmin1(k)=theoretical_Ynmin1(k)+Precursor_w_TMT_Matrix_for_mixing_ratio(i,j);
            end
        end
    end
end
%normalize matrix to 1 for positions that will be compared
theoretical_Ynmin1=theoretical_Ynmin1.*binary_abundand_Ynmin1_positions;
theoretical_Ynmin1=theoretical_Ynmin1./sum(theoretical_Ynmin1);

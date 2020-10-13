exp.well_quantified=zeros(size(exp.ratios,1),1);
for index=1:size(exp.ratios,1)
    if exp.ratios(index,3)>0.1
        exp.how_well_quantified(index) = sum(sqrt((exp.ratios(index,:)-exp.known_ratios).^2))/100;
    end
end

gscatter(exp.data.sum_ions_Ynmin1,exp.goodnes_fit(:), exp.how_well_quantified < 0.1);
set(gca,'XScale','log','Yscale','log','XLim',[100,max(exp.data.sum_ions_Ynmin1)],'YDir','reverse', 'YLim',[min(exp.goodnes_fit(:)),1])
legend('high quality: cosine distance < 0.05', 'poor quality: cosine distance > 0.05')
refline(0,exp.goodnesfit_cutoff)
y_range = get(gca,'Ylim');
line([exp.sum_ions_Ynmin1_cutoff, exp.sum_ions_Ynmin1_cutoff],[y_range(2),0.000000001]);
xlabel('Sum of ions in Pseudo Y(n-1) iso-envelope')
ylabel('Total square-difference of theoretical and measured Pseudo Y(n-1) iso-envelope')
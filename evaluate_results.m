% TGR_08996_TMTPro_Yeast-HeLa_8plex_Interference_1to1mix_ITrapid_z2_yeastonly.mat
load TGR_09069_TMTPro_Yeast-HeLa_8plex_Interference_1to1mix_ITrapid.mat
% TGR_08957_TMTPro_8plex_Lance_11_11_19_std.mat
% TGR_08959_TMTPro_8plex_Lance_11_11_19_wPatch.mat
% TGR_08961_TMTPro_8plex_Lance_11_11_19_wPatch_0_1shift.mat
% TGR_08801_HeLa_TMTPro8plex_ITMS_normal.mat
% TGR_08789_HeLa_TMTPro8plex_ITMS_enhanced.mat
% TGR_08787_HeLa_TMTPro_ITMS_zoom.mat
% TGR_08791_HeLa_TMTPro8plex_FTMS.mat
%plot_run_summaries(results)

data = read_in_data_csv('TGR_09069_TMTPro_Yeast-HeLa_8plex_Interference_1to1mix_ITrapid.xlsx');

%to_export = horzcat(results.ratios(:,2:9),results.data.sum_ions_Ynmin1, results.data.z);
%csvwrite('TGR_09069_TMTProcplus_output.csv',to_export)

%sum(sum(results.ratios<(median(results.ratios)/5)))
%sum(sum(results.ratios>(median(results.ratios)*5)))

% boxplot_percentwhisk_mod(results.ratios(results.indeces,:),95);

%data = read_in_data_csv('TGR_08957_TMTPro_8plex_Lance_11_11_19_std.xlsx');
%data = read_in_data_csv('TGR_08959_TMTPro_8plex_Lance_11_11_19_wPatch.xlsx');
%data = read_in_data_csv('TGR_08961_TMTPro_8plex_Lance_11_11_19_wPatch_0_1shift.xlsx');
%data = read_in_data_csv('TGR_08990_TMTPro_Yeast-HeLa_8plex_Interference_1to1mix_yeastonly.xlsx');

z3_cv = std(results.ratios(data.z==3 & results.indeces,2:9))./mean(results.ratios(data.z==3 & results.indeces,2:9));
z2_cv = std(results.ratios(data.z==2 & results.indeces,2:9))./mean(results.ratios(data.z==2 & results.indeces,2:9));



figure(1)
bar([z2_cv; z3_cv]')
    title('2+ and 3+ CV')
    xlabel('Channel')
    ylabel('CV')
    legend('2+','3+')

figure(2)
median(results.ratios(data.z==2 & results.data.sum_ions_Ynmin1>,2))



%sort(results.ratios(results.indeces,2))

figure(3)
histogram(results.ratios(data.z==2 & results.indeces, 9), 100)

%plot_run_summaries(results)
%{
figure
results.ratios(results.indeces,:);
how_many_MS_MS_quantified = sum(results.indeces & data.z==3);
how_many_unique_peps_quantified = size(unique((results.data.pep_sequences(results.indeces & data.z==3))),1);
how_many_protenis_quantified = size(unique((results.data.prot_name(results.indeces & data.z==3))),1);
how_many_original=size(results.ratios,1);
cleaned_identifier = strrep(results.name,'_',' ');
boxplot_percentwhisk_mod(results.ratios(results.indeces & data.z==3,:),95);
title_ = strcat(num2str(how_many_MS_MS_quantified),' of ',num2str(how_many_original),' MS/MS, ',num2str(how_many_unique_peps_quantified),' unique peps, ',num2str(how_many_protenis_quantified),' proteins');
title(title_)
ylabel('Relative abundance')
xtickangle(45);
set(gca,'Xtick',1:10,'Xticklabel',{'TMTPro0','TMTPro126','TMTPro127N','TMTPro129N','TMTPro130N','TMTPro131N','TMTPro131C','TMTPro132N','TMTPro134N','TMTPro135'})
%}

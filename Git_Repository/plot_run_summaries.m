%function to plot run summaries of Y(n-1) quantitation
function exp = plot_run_summaries(exp)

%Plot Boxplot of ratios
figure
subplot(3,3,1)
exp.ratios(exp.indeces,:);
how_many_MS_MS_quantified = sum(exp.indeces);
how_many_unique_peps_quantified = size(unique((exp.data.pep_sequences(exp.indeces))),1);
how_many_protenis_quantified = size(unique((exp.data.prot_name(exp.indeces))),1);
how_many_original=size(exp.ratios,1);
cleaned_identifier = strrep(exp.name,'_',' ');
boxplot_percentwhisk_mod(exp.ratios(exp.indeces,:),95);
title_ = strcat(num2str(how_many_MS_MS_quantified),' of ',num2str(how_many_original),' MS/MS, ',num2str(how_many_unique_peps_quantified),' unique peps, ',num2str(how_many_protenis_quantified),' proteins');
title(title_)
ylabel('Relative abundance')
xtickangle(45);
set(gca,'Xtick',1:10,'Xticklabel',{'TMTPro0','TMTPro126','TMTPro127N','TMTPro129N','TMTPro130N','TMTPro131N','TMTPro131C','TMTPro132N','TMTPro134N','TMTPro135'})

%Plot histograms for each channel
subplot(3,3,2)
n_all=[];
for index=1:size(exp.ratios,2)
    [n,xout]=hist(exp.ratios(exp.indeces,index),0:0.1:5);
    n_all=[n_all;n];
end
plot(xout,n_all)
legend('TMTPro0','TMTPro126','TMTPro127N','TMTPro129N','TMTPro130N','TMTPro131N','TMTPro131C','TMTPro132N','TMTPro134N','TMTPro135')
title_=strcat(cleaned_identifier,' Hist Ratios');
title(title_)

%Plot Histogram scatter plot of precursor intensity versus sum TMTc ions
subplot(3,3,3)
% Precursor intensity is dividec by 11.3 as Thermo seems to overestimat the
% number of ions (See Labnotebook entry from 2015-August-07)
loglog(exp.data.MS1_precursor_intensity(exp.data.mobile_protons<1).*exp.data.ion_injection_time(exp.data.mobile_protons<1)./1000./11.3,exp.data.sum_ions_Ynmin1(exp.data.mobile_protons<1),'.')
hold on
loglog(exp.data.MS1_precursor_intensity(exp.data.mobile_protons>=1).*exp.data.ion_injection_time(exp.data.mobile_protons>=1)./1000./11.3,exp.data.sum_ions_Ynmin1(exp.data.mobile_protons>=1),'rx')
xlabel('Injected Precursor (Intensity * Inj-time)')
ylabel('Sum charges in TMTc envelope (ions)')
legend({'no mobile protons','high mobility proton'})
% Calculate the median percentage of TMTc ions as a fraction of injected
fraction_TMTc_of_Precursor = 1/nanmedian(exp.data.MS1_precursor_intensity(exp.data.mobile_protons<1).*exp.data.ion_injection_time(exp.data.mobile_protons<1)./1000./11.3./exp.data.sum_ions_Ynmin1(exp.data.mobile_protons<1));
percent_conversion = num2str(fraction_TMTc_of_Precursor*100);
title({['Estim. ',percent_conversion,' % of precursor are converted into TMTc']})

%Plot Sum of Y(n-1)/charge for peptide with and without mobile proton
subplot(3,3,4)
n_all=[];
for index=0:1
    [n,xout]=hist(exp.data.sum_ions_Ynmin1(exp.mobile_protons_index == index),0:50:3000);
    n_all=[n_all;n];
    %plot(xout,n)
    %hold on
end
plot(xout,n_all)
legend('No mobile proton','Mobile proton')
%title_=strcat(cleaned_identifier,' Sum of ions in Y(n-1) envelope');
%title(title_)
y_range = get(gca,'Ylim');
line([exp.sum_ions_Ynmin1_cutoff, exp.sum_ions_Ynmin1_cutoff],[y_range(2), y_range(1)]);
xlabel('Number of ions in TMTc envelope')
ylabel('Number of peptides')



subplot(3,3,5)
%Plot a figure that explains cutoffs for goodness of fit

if exp.calculate_cosine_distance
    group = double(exp.cosine_distance < 0.01); % distinguish based on goodnesfit
    group(~exp.mobile_protons_index) = 3;        %mobile proton peptides get their own group
    gscatter(exp.data.sum_ions_Ynmin1,exp.goodnes_fit,group);
else
    group = zeros(size(exp.data.ScanF,1),1);
    group(~exp.mobile_protons_index) = 3;        %mobile proton peptides get their own group
    gscatter(exp.data.sum_ions_Ynmin1,exp.goodnes_fit,group);
end

%if exp.calculate_cosine_distance
%    goodnesfit_marker = exp.cosine_distance < 0.05;
%    gscatter(exp.data.sum_ions_Ynmin1(exp.indeces_yeast_only),exp.goodnes_fit(exp.indeces_yeast_only),goodnesfit_marker(exp.indeces_yeast_only));
%else
%    goodnesfit_marker = zeros(size(exp.data.ScanF,1),1);
%    gscatter(exp.data.sum_ions_Ynmin1,exp.goodnes_fit(:),goodnesfit_marker);
%end


set(gca,'XScale','log','Yscale','log','XLim',[50,max(exp.data.sum_ions_Ynmin1)],'YDir','reverse', 'YLim',[min(exp.goodnes_fit(:)),1])
legend('Mobile Proton','poor quality: cosine distance > 0.05','high quality: cosine distance < 0.05')
refline(0,exp.goodnesfit_cutoff)
%plot function for quality cutof
hold on
if isfield(exp, 'a_quality_cutoff') && isfield(exp, 'b_quality_cutoff')
   fnch = @(x) exp.a_quality_cutoff*(x^exp.b_quality_cutoff);
    fplot(fnch,[50,max(exp.data.sum_ions_Ynmin1)])    
end

%plot line indicating Sum TMTc cluster cutoff
y_range = get(gca,'Ylim');
line([exp.sum_ions_Ynmin1_cutoff, exp.sum_ions_Ynmin1_cutoff],[y_range(2),0.000000001]);
xlabel('Sum of ions in TMTc iso-envelope')
ylabel('Total square-difference of theoretical and measured Pseudo Y(n-1) iso-envelope')

%Plot Sum of Y(n-1)/charge for each charge state
subplot(3,3,6)
n_all=[];
for index=2:max(exp.data.z)
    [n,xout]=hist(exp.data.sum_ions_Ynmin1(exp.data.z==index),0:10:4000);
    n_all=[n_all;n];
    %plot(xout,n)
    %hold on
end
plot(xout,n_all)
legend('2+','3+','4+','5+','6+')
title_=strcat(cleaned_identifier,' Sum of ions in Y(n-1) envelope');
title(title_)
y_range = get(gca,'Ylim');
line([exp.sum_ions_Ynmin1_cutoff, exp.sum_ions_Ynmin1_cutoff],[y_range(2), y_range(1)]);




%Plot correlation bw Precursor intensity (in spectrum) and Ynmin1 intensity
subplot(3,3,7)
if any(strcmp('precursor_iso_envelope',fieldnames(exp.data)))
    exp.data.sum_precursor_envelope = sum(exp.data.precursor_iso_envelope,2);
    plot(exp.data.sum_ions_Ynmin1(~exp.mobile_protons_index & exp.data.z == 2),exp.data.sum_precursor_envelope(~exp.mobile_protons_index & exp.data.z == 2),'bx')
    hold on
    plot(exp.data.sum_ions_Ynmin1(~exp.mobile_protons_index & exp.data.z == 3),exp.data.sum_precursor_envelope(~exp.mobile_protons_index & exp.data.z == 3),'rx')
    set(gca,'Xlim',[0,4000],'Ylim',[0,4000]);
    lsline
    xlabel('Sum TMTc ions')
    ylabel('Sum Surving precursor ions')
    title('No mobile proton')
    legend({'2+','3+'})
end


%Plot correlation bw summ of MS2 low m/z reporter ions and summed TMTc
%intensity
subplot(3,3,8)
plot(exp.data.sum_ions_Ynmin1(~exp.mobile_protons_index),sum(exp.data.reporter_ions(~exp.mobile_protons_index,:),2),'x')
set(gca,'Xlim',[0,4000],'Ylim',[0,15000]);
xlabel('Sum TMTc ions')
ylabel('Sum Reporter ions')
title('No mobile proton')


if exp.calculate_cosine_distance
    %Calculate FDR
    exp.num_quantified=sum(exp.indeces);
    exp.Quant_false=sum(exp.indeces & exp.cosine_distance > 0.02);
    exp.Quant_FDR = exp.Quant_false/sum(exp.indeces);
end

if  any(strcmp('MS1_precursor_intensity',fieldnames(exp.data)))
    subplot(3,3,9)
    [n,xout] = hist(log10(exp.data.MS1_precursor_intensity(exp.indeces & exp.data.z ==2)),[3:0.1:10]);
    plot(xout,n,'x-r')
    hold on
    [n,xout] = hist(log10(exp.data.MS1_precursor_intensity(exp.data.z ==2)),[3:0.1:10]);
    plot(xout,n,'x-b')
    [n,xout] = hist(log10(exp.data.MS1_precursor_intensity(exp.indeces & exp.data.z ==3)),[3:0.1:10]);
    plot(xout,n,'x-g')
    hold on
    [n,xout] = hist(log10(exp.data.MS1_precursor_intensity(exp.data.z ==3)),[3:0.1:10]);
    plot(xout,n,'x-k')
    legend({'Succesfully quantified 2+','All 2+','Succesfully quantified 3+','All 3+'})
    xlabel('log10 Precursor Intensity')
end

end

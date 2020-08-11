%Calcualte and plot ratios for run 
function exp=plot_ratios_goodnes_fit(exp)
%Define binsize for variance plot
binsize=0.0001; %binsize for SumYions
%Avoid division by zero
exp.ratios_of_channels = zeros(size(exp.ratios,1),5);
%make colormap
colormap = hsv(5);
%Calculate max bin
max_bin = round(max(exp.goodnes_fit)/binsize);

for index_column=[1,2,4,5]
    %go throug structure and calculate channels ratios
    for index_row = 1:size(exp.ratios,1)
        %Avoid division by 0
        if exp.ratios(index_row,3)> 0.01;
           exp.ratios_of_channels(index_row,index_column) = exp.ratios(index_row,index_column)/exp.ratios(index_row,3);
        end
    end
    %plot(exp.ratios_of_channels(:,index_column), exp.data.sum_ions_Ynmin1,'x','Color',colormap(index_column,:))
    %set(gca,'Color',colormap(index_column,:));
    %hold on 
    %set(gca,'Xlim',[0,2]);

    % Calculate known_channel_ratios
    exp.known_ratios_of_channels(index_column) = exp.known_ratios(index_column)/exp.known_ratios(3);

    for n = 1:max_bin
            index_for_bin = exp.goodnes_fit > (n-1).*binsize & exp.goodnes_fit <= n.*binsize & exp.data.sum_ions_Ynmin1 >1;
            median_variance(n) = median(sqrt((exp.ratios_of_channels(index_for_bin,index_column)-exp.known_ratios_of_channels(index_column)).^2));
    end
    plot((0.5:max_bin).*binsize,median_variance,'x','Color',colormap(index_column,:))
    hold on 
end
legend('1:10 w interference','4:10 w interference','4:10 wo interference','1:10 wo interference')
set(gca,'Ylim',[0,0.5],'Xlim',[0,0.01]);
xlabel('summed squared Diff theoretical and observed Y(n-1) envelope');
ylabel('Median Variance');
end



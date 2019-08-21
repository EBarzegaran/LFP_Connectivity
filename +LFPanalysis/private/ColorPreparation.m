% this script prepare the colors to match the paper: Plomp et al., EJN 2014
% requires brewermap toolbox
clear; clc;close all
Colors = brewermap(6,'Set2'); % Set colors for layers 
Colors = Colors([2 6 5 1 4 3],:);

for c = 1:6
    plot(rand(30,1)+.1*c,'Color',Colors(c,:),'linewidth',1.5)
    hold on;
end

legend(arrayfun(@(x) ['L' num2str(x)],1:6,'uni',false));

save('LayerColors.mat','Colors');


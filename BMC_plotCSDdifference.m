clear

close all
cd 'E:\LaCie\DATA_KD\BRFS_fCSD'
load fCSDfrom161007_001.mat
Session = '161007_E_001';


difference = LargeResp - SmallResp;
Figuresdir = 'E:\LaCie\DATA_KD\BRFS_figs';
%% Plot CSD for stim 1, no soa, dCOS - bi
pre = 100;
post = 500;
tvec = (-pre:post);
chanvec = linspace(1,24,size(difference,1));



figure

imagesc(tvec,chanvec,difference); 
colormap(gray);  colorbar; vline(0);
set(gca,'XTickLabel',[700 800 900 1000 1100 1200 1300]);
mn2 = min(min(difference)); mx2 = max(max(difference)); 
maxval2 = max([abs(mn2) abs(mx2)]);
 caxis([-maxval2 maxval2]);
% % % % caxis([-1000 1000]);
title('161007 E 001 difference. 75-165. Sink Btm 17');
xlabel('time (ms)') 
ylabel('electrode number') 

saveas(gcf,fullfile(Figuresdir,strcat(Session,'_difference')), 'jpeg')


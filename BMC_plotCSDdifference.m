
close all
cd 'E:\LaCie\DATA_KD\161005_E\'



%% Plot CSD
pre = 100;
post = 500;
tvec = (-pre:post);
chanvec_elC = linspace(2,23,size(diff_elC,2));
chanvec_elD = linspace(2,24,size(diff_elC,2));

figure(3)
subplot(1,2,1); 
imagesc(tvec,chanvec_elC,diff_elC);
colormap(flipud(jet)); colorbar;  vline(0);
set(gca,'XTickLabel',[800 1000 1200]);
mn1 = min(min(diff_elC)); mx1 = max(max(diff_elC));
maxval1 = max([abs(mn1) abs(mx1)]);
caxis([-maxval1 maxval1]);
title('161005 el 1 dCOS-Bi soa stim 2');
xlabel('time (ms)') 
ylabel('electrode number') 

subplot(1,2,2);
imagesc(tvec,chanvec_elD,diff_elD); 
colormap(flipud(jet));  colorbar; vline(0);
set(gca,'XTickLabel',[800 1000 1200]);
mn2 = min(min(diff_elD)); mx2 = max(max(diff_elD)); 
maxval2 = max([abs(mn2) abs(mx2)]);
% % % % caxis([-maxval2 maxval2]);
caxis([-300 300]);
title('161005 el 2 dCOS-bi soa stim2');
xlabel('time (ms)') 
ylabel('electrode number') 

% % % % dateFormatOut = 'yyyy-mm-dd';
% % % % saveDate = datestr(now,dateFormatOut);
% % % % saveName = strcat('fig_CSD_',Filename,'_',saveDate);
% % % % clear dateFormatOut saveDate
% % % % saveas(figure(3),saveName)


%CONGRATS... it finished
load gong
sound(y,Fs)


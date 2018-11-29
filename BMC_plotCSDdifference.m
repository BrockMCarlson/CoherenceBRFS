clear

close all
cd 'E:\LaCie\DATA_KD\161005_E\'

load BRFS_bi_1.mat
load BRFS_dCOS_1.mat
load BRFS_bi_2.mat
load BRFS_dCOS_2.mat

stim1.elC = dCOS_1_elC - bi_1_elC;
stim1.elD = dCOS_1_elD - bi_1_elD;

stim2.elC = dCOS_2_elC - bi_2_elC;
stim2.elD = dCOS_2_elD - bi_2_elD;

clearvars -except stim1 stim2
%% Plot CSD for stim 1, no soa, dCOS - bi
pre = 100;
post = 500;
tvec = (-pre:post);
chanvec_1_elC = linspace(2,23,size(stim1.elC,2));
chanvec_1_elD = linspace(2,24,size(stim1.elD,2));

% elC
figure(1)
subplot(1,2,1); 
imagesc(tvec,chanvec_1_elC,stim1.elC);
colormap(flipud(jet)); colorbar;  vline(0);
set(gca,'XTickLabel',[0 200 400]);
mn1 = min(min(stim1.elC)); mx1 = max(max(stim1.elC));
maxval1 = max([abs(mn1) abs(mx1)]);
caxis([-maxval1 maxval1]);
title('161005 el 1 dCOS-Bi stim 1');
xlabel('time (ms)') 
ylabel('electrode number') 

%elD
subplot(1,2,2);
imagesc(tvec,chanvec_1_elD,stim2.elD); 
colormap(flipud(jet));  colorbar; vline(0);
set(gca,'XTickLabel',[0 200 400]);
mn2 = min(min(stim2.elD)); mx2 = max(max(stim2.elD)); 
maxval2 = max([abs(mn2) abs(mx2)]);
% % % % caxis([-maxval2 maxval2]);
caxis([-100 100]);
title('161005 el 2 dCOS-bi stim 1');
xlabel('time (ms)') 
ylabel('electrode number') 


%% Plot CSD for stim 2, after soa, dCOS - bi
pre = 100;
post = 500;
tvec = (-pre:post);
chanvec_2_elC = linspace(2,23,size(stim2.elC,2));
chanvec_2_elD = linspace(2,24,size(stim2.elD,2));


figure(2)
subplot(1,2,1); 
imagesc(tvec,chanvec_2_elC,stim2.elC);
colormap(flipud(jet)); colorbar;  vline(0);
set(gca,'XTickLabel',[800 1000 1200]);
mn1 = min(min(stim2.elC)); mx1 = max(max(stim2.elC));
maxval1 = max([abs(mn1) abs(mx1)]);
caxis([-maxval1 maxval1]);
title('161005 el 1 dCOS-Bi soa stim 2');
xlabel('time (ms)') 
ylabel('electrode number') 

subplot(1,2,2);
imagesc(tvec,chanvec_2_elD,stim2.elD); 
colormap(flipud(jet));  colorbar; vline(0);
set(gca,'XTickLabel',[800 1000 1200]);
mn2 = min(min(stim2.elD)); mx2 = max(max(stim2.elD)); 
maxval2 = max([abs(mn2) abs(mx2)]);
caxis([-maxval2 maxval2]);
% % % % caxis([-200 200]);
title('161005 el 2 dCOS-bi soa stim2');
xlabel('time (ms)') 
ylabel('electrode number') 


%% CONGRATS... it finished
load gong
sound(y,Fs)


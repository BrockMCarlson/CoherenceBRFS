%BMC plot CSD
%loads a matlab workspace created in ReadInNS6 to use variables for
%processing and then plotting. Available variables, 'EV.tp' 'Filename'
%'LFP' 'MUA' and saveName'
%trigger to stimulus on based on event codes
%average LFP across trials
    %run calcCSD
    % convert uV to nA/mm^3
    % subtract baseline
    % 1/2 wave rectify
    % pad out to same number of channels as LFP
    % filter CSD
    %plot

clear
close all

cd 'E:\LaCie\DATA_KD\brfs\Brock\151222_E'
%load WorkspaceForProcessing_180705_evp001_test.mat
load BRFSWorkspace001_151222_E_brfs001ns2_2018-11-26.mat
% Available variables, 'Cond' 'EV' 'extension' 'Filename' 'LFP' 'saveName'
% and 'Unq_cond' 


%% 4. PROCESS
clear stimLFP
pre = 100;   % pre-stim time (baseline) in ms
post = 500; % post-stim time in ms

% TRIGGER MUA & LFP TO STIM ON
for tr = 1:length(EV.D) % trigger to stim-on times for all trials
    stimtm = round(EV.D(tr,2)/30) ;% divide by 30 to convert to 1kHz. Note, LFP and MUA already in 1kHZ
    refwin = stimtm-pre:stimtm+post;
    stimLFP(tr,:,:) = LFP(refwin,:);
end

% COMPUTE AVERAGE ACROSS TRIALS
avgLFP = squeeze(mean(stimLFP,1));
avgLFP_elC = avgLFP(:,1:24);
avgLFP_elD = avgLFP(:,25:48);

tvec = (-pre:post);
figure(1),cla
title('161005 elC CSDline');
    for chanC = 1:size(avgLFP_elC,2)
       plC = subplot(24,1,chanC);plot(tvec,avgLFP_elC(:,chanC));vline(0);hline(0);  
    end
    

    figure(2),
    title('161005 elD CSDline');
    for chanD = 1:size(avgLFP_elD,2)
        subplot(24,1,chanD);plot(tvec,avgLFP_elD(:,chanD));vline(0);hline(0);
    end

%% Get CSD
CSD_elC = calcCSD(avgLFP_elC).*0.4; 
CSD_elD = calcCSD(avgLFP_elD).*0.4;

bl =1:50;
bl_CSD_elC = mean(CSD_elC(:,bl),2);
bl_CSD_elD = mean(CSD_elD(:,bl),2);
blCSD_elC = (CSD_elC-bl_CSD_elC); 
blCSD_elD = (CSD_elD-bl_CSD_elD);

gauss_sigma = 0.1;
padCSD_elC = padarray(blCSD_elC,[1 0],NaN,'replicate');
padCSD_elD = padarray(blCSD_elD,[1 0],NaN,'replicate');
fCSD_elC = filterCSD(padCSD_elC,gauss_sigma);
fCSD_elD = filterCSD(padCSD_elD,gauss_sigma);

%% Plot CSD
tvec = (-pre:post);
chanvec_elC = linspace(2,23,size(fCSD_elC,2));
chanvec_elD = linspace(2,24,size(fCSD_elD,2));

figure(3)
subplot(1,2,1); 
imagesc(tvec,chanvec_elC,fCSD_elC);
colormap(flipud(jet)); colorbar;  vline(0);
set(gca,'XTickLabel',[800 1000 1200]);
mn1 = min(min(fCSD_elC)); mx1 = max(max(fCSD_elC));
maxval1 = max([abs(mn1) abs(mx1)]);
caxis([-maxval1 maxval1]);
title('161005 el 1 dCOS soa stim2');
xlabel('time (ms)') 
ylabel('electrode number') 

subplot(1,2,2);
imagesc(tvec,chanvec_elD,fCSD_elD); 
colormap(flipud(jet));  colorbar; vline(0);
set(gca,'XTickLabel',[800 1000 1200]);
mn2 = min(min(fCSD_elD)); mx2 = max(max(fCSD_elD)); 
maxval2 = max([abs(mn2) abs(mx2)]);
caxis([-maxval2 maxval2]);
% % % %  caxis([-300 300]);
title('161005 el 2 dCOS soa stim2');
xlabel('time (ms)') 
ylabel('electrode number') 

% % % % dateFormatOut = 'yyyy-mm-dd';
% % % % saveDate = datestr(now,dateFormatOut);
% % % % saveName = strcat('fig_CSD_',Filename,'_',saveDate);
% % % % clear dateFormatOut saveDate
% % % % saveas(figure(3),saveName)

dCOS_2_elC = fCSD_elC;
dCOS_2_elD = fCSD_elD;
save BRFS_dCOS_2.mat dCOS_2_elC dCOS_2_elD

%CONGRATS... it finished
load gong
sound(y,Fs)


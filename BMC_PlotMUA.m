%BMC getStarted 180627 Process and Plot MUA&LFP
%loads a matlab workspace created in MakeWorspace to use variables for
%processing and then plotting
%trigger to stimulus on based on event codes
%average across trials
%plot
clear
close all
cd 'E:\LaCie\DATA_KD\160831_E\'

%load WorkspaceForProcessing_180705_evp001_test.mat
load Workspace_160831_E_evp002_2018-08-02.mat
% Available variables, 'EV.tp' 'Filename' 'LFP' 'MUA' and saveName'

%% 4. PROCESS
clear EVP trigLFP
pre = 100;   % pre-stim time (baseline) in ms
post = 500; % post-stim time in ms

% TRIGGER MUA & LFP TO STIM ON
for tr = 1:length(EV.tp) % trigger MUA to stim-on times for all trials
    stimtm = round(EV.tp(tr,1)/30) ; % divide by 30 to convert to 1kHz
    refwin = stimtm-pre:stimtm+post;
    EVP(tr,:,:) = MUA(refwin,:);
    stimLFP(tr,:,:) = LFP(refwin,:);
end

 

% COMPUTE AVERAGE ACROSS TRIALS
xVEP = squeeze(mean(EVP,1));
VEP = xVEP;   %(:,ids);
avgLFP = squeeze(mean(stimLFP,1));
%avgLFP = xLFP;  %(:,ids)   %./4; % divide by 4 here?? To get to mV


%% 5. PLOT
tvec = (-pre:post);
figure(1),cla
    for chan = 1:size(VEP,2)
        sp1_1 = subplot(1,2,1);plot(tvec,VEP(:,chan)-0.5*chan);vline(0);
        hold on
    end
axis tight

    for chan = 1:size(avgLFP,2)
        sp1_2 = subplot(1,2,2);plot(tvec,avgLFP(:,chan)-0.5*chan) ;vline(0);
        hold on
    end
axis tight

% baseline correct
bl =pre;
bl_mua = mean(VEP(bl,:),1);
blVEP = (VEP-bl_mua);    %./bl_mua;
bl_LFP = mean(avgLFP(bl,:),1);
blLFP = (avgLFP-bl_LFP);    %./bl_LFP;

% plot in 2D
figure(2)
tvec = (-pre:post);
sp2_1 = subplot(2,2,1); imagesc(tvec,[],blVEP(:,1:24)'); colorbar; vline(0);
mn1 = min(min(blVEP(:,1:24))) .* .7; mx1 = max(max(blVEP(:,1:24))) .* .7; 
caxis([mn1 mx1]); 
sp2_2 = subplot(2,2,2); imagesc(tvec,[],blVEP(:,25:48)'); colorbar; vline(0);
mn2 = min(min(blVEP(:,25:48))) .* .7; mx2 = max(max(blVEP(:,25:48))) .* .7; 
caxis([mn2 mx2]); 
sp2_3 = subplot(2,2,3); imagesc(tvec,[],blLFP(:,1:24)'); colorbar; vline(0);
mn3 = min(min(blLFP(:,1:24))) .* .7; mx3 = max(max(blLFP(:,1:24))) .* .7; 
caxis([mn3 mx3]); 
sp2_4 = subplot(2,2,4); imagesc(tvec,[],blLFP(:,25:48)'); colorbar; vline(0);
mn4 = min(min(blLFP(:,25:48))) .* .7; mx4 = max(max(blLFP(:,25:48))) .* .7; 
caxis([mn4 mx4]); 

title(sp1_1,'VEP across channels');
title(sp1_2,'LFP across channels');
title(sp2_1,'% Change, VEP, electrode 1');
title(sp2_2,'% Change, VEP, electrode 2');
title(sp2_3,'% Change, LFP, electrode 1');
title(sp2_4,'% Change, LFP, electrode 2');

dateFormatOut = 'yyyy-mm-dd';
saveDate = datestr(now,dateFormatOut);
saveName = strcat('fig_MUA_',Filename,'_',saveDate);
clear dateFormatOut saveDate



saveas(figure(1),strcat(saveName,'_1'))
saveas(figure(2),strcat(saveName,'_2'))

%CONGRATS... it finished
load gong
sound(y,Fs)

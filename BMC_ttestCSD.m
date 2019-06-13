%BMC CSD ttest
clear
close all
cd 'E:\LaCie\DATA_KD\161005_E\'
%load WorkspaceForProcessing_180705_evp001_test.mat
load BRFSWorkspace_161005_E_brfs001ns2_2018-11-12.mat
% Available variables, 'Cond' 'EV' 'extension' 'Filename' 'LFP' 'saveName'
% and 'Unq_cond' 

%% 4. PROCESS
clear stimLFP
pre = 100;   % pre-stim time (baseline) in ms
post = 500; % post-stim time in ms

% TRIGGER MUA & LFP TO STIM ON
for tr = 1:length(EV.B) % trigger to stim-on times for all trials
    stimtm = round(EV.B(tr,1)/30) ;% divide by 30 to convert to 1kHz. Note, LFP and MUA already in 1kHZ
    refwin = stimtm-pre:stimtm+post;
    stimLFP(tr,:,:) = LFP(refwin,:);
end

switchOrderStimLFP = permute(stimLFP,[2 3 1]);
LFP_elC = switchOrderStimLFP(:,1:24,:);
LFP_elD = switchOrderStimLFP(:,25:48,:);

%% Get CSD
CSD_elC = calcCSD(LFP_elC).*0.4; 
CSD_elD = calcCSD(LFP_elD).*0.4;

gauss_sigma = 0.1;
padCSD_elC = padarray(CSD_elC,[1 0],NaN,'replicate');
padCSD_elD = padarray(CSD_elD,[1 0],NaN,'replicate');

%% t-test each point in time after stimulus onset on each electrode, t-test
% accross the trials.

clear i j elC elD h
collTstat = zeros(601,1);
tResult = zeros(24,601);
for i = 1:size(A.C,1) % performed per electrode contact
    
    contact_Bi = squeeze(A.D(i,:,:));
    contact_dCOS = squeeze(B.D(i,:,61:122));
    
    for j = 1:size(contact_Bi,1)
        Bi = contact_Bi(j,:);
        dCOS = contact_dCOS(j,:);
        [h,p,ci,stats] = ttest(Bi,dCOS);
        collTstat(j) = stats.tstat;  
    
    end
    
    tResult(i,:) =  collTstat;
    
end

%% Plot tResult
pre = 100;   % pre-stim time (baseline) in ms
post = 500; % post-stim time in ms
tvec = (-pre:post);
chanvec = size(tResult,1);

figure
imagesc(tvec,1,tResult);
colormap(jet); colorbar; vline(0);
mn1 = min(min(tResult)); mx1 = max(max(tResult));
maxval1 = max([abs(mn1) abs(mx1)]);
caxis([-maxval1 maxval1]);
title('t scores dCOS vs Bi on electrode D. 2nd half dCOS sample');


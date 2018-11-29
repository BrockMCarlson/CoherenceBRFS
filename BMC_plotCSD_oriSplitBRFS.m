% BMC_plotCSD_oriSplitBRFS
% November 27th, 2018.
clear
close all
% This script creates 3 figures with 2 subplots each. Together they will 
% form a figure with the first row focused on 75deg presentation. 
% The second row is focused on 165deg presentation. 

% Only conditions where matching contrast of .9 in both eyes is examined.

% In the first column (fig 1; 2,1,1 & 2,1,2) are the CSD line plots for one
% channel above the sink bottom. These are found for binocular (dioptic)
% simultaneous (no soa) presentations of parallel gratings at 75deg and
% 165deg, respectively. This comes from condition group A.
% 75  controll -- Cond.A = [: : 75  75  0.9 0.9 0 2];
% 165 controll -- Cond.A = [: : 165 165 0.9 0.9 0 2];

% In the first row, second and third column (fig 2; 1,2,1 & ,1,2,2), Are the CSD
% plots for condition group D when 75 degrees orientation is presented
% on the flash. 75deg grating may appear in either eye, as long as it is presented second (as the flash).
% The subplot 1 is the filtered CSD across all channels. The subplot 2 is
% the line plot of channels a)5-above-sink b)sink c)5-below-sink

% % In the second row, second and third column (fig 3; 1,2,1 & ,1,2,2), Are the
% CSD plots for condiion grou D when 165 deg ori is presented on the flash.
% Subplot 4 = filtered CSD across all channels. Subplot 5 is line plot of
% channel a)5-above-sink b)sink and c)5-below-sink

% Needed input: Cond, EV, and LFP
% Change every time 
%   sink
%   titles

Sink = 15;
%% 1. Load data

cd 'E:\LaCie\DATA_KD\161005_E\'

%load WorkspaceForProcessing_180705_evp001_test.mat
load BRFSWorkspace_161005_E_brfs001ns2_2018-11-16.mat

%% 2. Plot Controll - Simultaneous Binocular Dioptic Grating 

%Create EV.A subsets - "biOri75" and "biOri165"
clear count75 count165
a = 0;
count165 = 0;
for i = 1:length(Cond.A)
    if Cond.A(i,3:8) == [75,75,0.9,0.9,0,2]
        a = a+1;
        biOri75(a,:) = EV.A(i,:);
    elseif Cond.A(i,3:8) == [165,165,0.9,0.9,0,2]
        count165 = count165+1;
        biOri165(count165,:) = EV.A(i,:);
    end
end

% Calculate CSD for controll
clear stimLFP
pre = 100;
post = 500;

% Trigger LFP to stim on
for j = 1:length(biOri75)
    stimtm75 = round(biOri75(j,1)/30); %divide by 30 to convert timepoints to 1kHz. LFP already in 1kHz
    refwin75 = stimtm75-pre:stimtm75+post;
    stimLFP75(j,:,:) = LFP(refwin75,:); 
end

for k = 1:length(biOri165)
    stimtm165 = round(biOri165(k,1)/30); %divide by 30 to convert timepoints to 1kHz. LFP already in 1kHz
    refwin165 = stimtm165-pre:stimtm165+post;
    stimLFP165(k,:,:) = LFP(refwin165,:); 
end

% Compute average across trials
avgLFP75 = squeeze(mean(stimLFP75,1));
avgLFP165 = squeeze(mean(stimLFP165,1));

% Calculate CSD and subtract baseline
CSD75 = calcCSD(avgLFP75(:,1:24)).*0.4;
CSD165 = calcCSD(avgLFP165(:,1:24)).*0.4;
bl =1:50;
bl_CSD75 = mean(CSD75(:,bl),2);
bl_CSD165 = mean(CSD165(:,bl),2);
blCSD75 = (CSD75-bl_CSD75); 
blCSD165 = (CSD165-bl_CSD165);

%% FIX

 %%%%  PLOT THE STANDARD ERROR ON THE MEAN ?

%% 

% Plot controll
chan = Sink-1; %this probably actually is 2 above the sink because when the... 
             %CSD is calculated, a channel is lost on top and bottom. 
tvec = -pre:post;
figure(1),cla
	h1= plot(tvec,blCSD75(chan,:)); hold on
	h2 = plot(tvec,blCSD165(chan,:));
        vline(0); 
        title('161005 ori 75vs165 CSD sink lines');        
        ylabel( 'nA/(mm^3)') ; 
        xlabel('time (ms)');
        legend([h1 h2],{'75','165'},'location','best');
        hold off
clearvars -except LFP EV Cond Sink
        
%% 3. Figure set 2, 75deg on flash, (second presentation) CSD

pre = 100;
post = 500;

%Create EV.D subset - dCOS75
a = 0;
for i = 1:length(Cond.D)
    if Cond.D(i,3:8) == [165,75,0.9,0.9,800,3]
        a = a+1;
        dCOSori75(a,:) = EV.D(i,:);
    end
end

% Trigger LFP to stim on
for j = 1:length(dCOSori75)
    stimtm_dCOSori75 = round(dCOSori75(j,1)/30); %divide by 30 to convert timepoints to 1kHz. LFP already in 1kHz
    refwin_dCOSori75 = stimtm_dCOSori75-pre:stimtm_dCOSori75+post;
    stimLFP_dCOSori75(j,:,:) = LFP(refwin_dCOSori75,:); 
end
 
% Compute average across trials
avgLFP_dCOSori75 = squeeze(mean(stimLFP_dCOSori75,1));

% Calculate CSD and subtract baseline
CSD_dCOSori75 = calcCSD(avgLFP_dCOSori75(:,1:24)).*0.4;
bl =1:50;
bl_CSD_dCOSori75 = mean(CSD_dCOSori75(:,bl),2);
blCSD_dCOSori75 = (CSD_dCOSori75-bl_CSD_dCOSori75); 

% pad 
padCSD_dCOSori75 = padarray(blCSD_dCOSori75,[1 0],NaN,'replicate'); %48x1301

% pull out line plots for later
    %supragranular 
    supGranChan = Sink - 7;
    supGran_75 = padCSD_dCOSori75(supGranChan,:);
    
    %granular
    granChan = Sink - 2;
    Gran_75 = padCSD_dCOSori75(granChan,:);
    
    %infragranular
    infraGranChan = Sink + 5;
    infraGran_75 = padCSD_dCOSori75(infraGranChan,:);
    
% filter
gauss_sigma = 0.1;
fCSD_dCOSori75 = filterCSD(padCSD_dCOSori75,gauss_sigma);

% plot filtered data
tvec = (-pre:post);
chanvec_dCOSori75 = linspace(2,24,size(fCSD_dCOSori75,1));

figure(2)

    imagesc(tvec,chanvec_dCOSori75,fCSD_dCOSori75);
        colormap(flipud(jet)); colorbar;  vline(0);
        set(gca,'XTickLabel',[700 800 900 1000 1100 1200 1300]);
        title('161005 dCOS 75 deg orientation');
        xlabel('time (ms)') 
        ylabel('electrode number') 
        mn1 = min(min(fCSD_dCOSori75)); mx1 = max(max(fCSD_dCOSori75));
        maxval1 = max([abs(mn1) abs(mx1)]);
        caxis([-maxval1 maxval1]);

figure(3)        
    ax1 = subplot(3,1,1);
    	plot(supGran_75);
        ylim([-2000 500])
        title('dCOS 75 supragranular');
    ax2 = subplot(3,1,2);
        plot(Gran_75);
        ylim([-2000 500])
        title('dCOS 75 Granular');
    ax3 = subplot(3,1,3);
        plot(infraGran_75)
        ylim([-2000 500])
        title('dCOS 75 infragranular');

        clearvars -except LFP EV Cond Sink
        
%% 4. Figure set 3, 165deg on flash, (second presentation) CSD

pre = 100;
post = 500;

%Create EV.D subset - dCOS75
a = 0;
for i = 1:length(Cond.D)
    if Cond.D(i,3:8) == [75,165,0.9,0.9,800,3]
        a = a+1;
        dCOSori165(a,:) = EV.D(i,:);
    end
end

% Trigger LFP to stim on
for j = 1:length(dCOSori165)
    stimtm_dCOSori165 = round(dCOSori165(j,1)/30); %divide by 30 to convert timepoints to 1kHz. LFP already in 1kHz
    refwin_dCOSori165 = stimtm_dCOSori165-pre:stimtm_dCOSori165+post;
    stimLFP_dCOSori165(j,:,:) = LFP(refwin_dCOSori165,:); 
end
 
% Compute average across trials
avgLFP_dCOSori165 = squeeze(mean(stimLFP_dCOSori165,1));

% Calculate CSD and subtract baseline
CSD_dCOSori165 = calcCSD(avgLFP_dCOSori165(:,1:24)).*0.4;
bl =1:50;
bl_CSD_dCOSori165 = mean(CSD_dCOSori165(:,bl),2);
blCSD_dCOSori165 = (CSD_dCOSori165-bl_CSD_dCOSori165); 

% pad 
padCSD_dCOSori165 = padarray(blCSD_dCOSori165,[1 0],NaN,'replicate'); %48x1301

% pull out line plots for later
    %supragranular 
    supGranChan = Sink - 7;
    supGran_75 = padCSD_dCOSori165(supGranChan,:);
    
    %granular
    granChan = Sink - 2;
    Gran_75 = padCSD_dCOSori165(granChan,:);
    
    %infragranular
    infraGranChan = Sink + 5;
    infraGran_75 = padCSD_dCOSori165(infraGranChan,:);
    
% filter
gauss_sigma = 0.1;
fCSD_dCOSori165 = filterCSD(padCSD_dCOSori165,gauss_sigma);

% plot filtered data
tvec = (-pre:post);
chanvec_dCOSori165 = linspace(2,24,size(fCSD_dCOSori165,1));

figure(4)

    imagesc(tvec,chanvec_dCOSori165,fCSD_dCOSori165);
        colormap(flipud(jet)); colorbar;  vline(0);
        set(gca,'XTickLabel',[700 800 900 1000 1100 1200 1300]);
        title('161005 dCOS 165 deg orientation');
        xlabel('time (ms)') 
        ylabel('electrode number') 
        mn1 = min(min(fCSD_dCOSori165)); mx1 = max(max(fCSD_dCOSori165));
        maxval1 = max([abs(mn1) abs(mx1)]);
        caxis([-maxval1 maxval1]);

figure(5)        
    ax1 = subplot(3,1,1);
    	plot(supGran_75);
        ylim([-2000 500])
        title('dCOS 165 supragranular');
    ax2 = subplot(3,1,2);
        plot(Gran_75);
        ylim([-2000 500])
        title('dCOS 165 Granular');
    ax3 = subplot(3,1,3);
        plot(infraGran_75)
        ylim([-2000 500])
        title('dCOS 165 infragranular');

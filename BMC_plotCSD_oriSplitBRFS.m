% BMC_plotCSD_oriSplitBRFS
% November 27th, 2018.
clear
close all
% This script creates 3 figures with 2 subplots each. Together they will 
% form a figure with the first row focused on Xdeg presentation. 
% The second row is focused on Ydeg presentation. 

% Only conditions where matching contrast of .9 in both eyes is examined.

% In the first column (fig 1; 2,1,1 & 2,1,2) are the CSD line plots for one
% channel above the sink bottom. These are found for binocular (dioptic)
% simultaneous (no soa) presentations of parallel gratings at Xdeg and
% Ydeg, respectively. This comes from condition group A.
% X  controll -- Cond.A = [: : X  X  0.9 0.9 0 2];
% Y controll -- Cond.A = [: : Y Y 0.9 0.9 0 2];

% In the first row, second and third column (fig 2; 1,2,1 & ,1,2,2), Are the CSD
% plots for condition group D when X degrees orientation is presented
% on the flash. Xdeg grating may appear in either eye, as long as it is presented second (as the flash).
% The subplot 1 is the filtered CSD across all channels. The subplot 2 is
% the line plot of channels a)5-above-sink b)sink c)5-below-sink

% % In the second row, second and third column (fig 3; 1,2,1 & ,1,2,2), Are the
% CSD plots for condiion grou D when Y deg ori is presented on the flash.
% Subplot 4 = filtered CSD across all channels. Subplot 5 is line plot of
% channel a)5-above-sink b)sink and c)5-below-sink

% Needed input: Cond, EV, and LFP
% Change every time 
%   sink
%   titles

Sink = 17;
Session = '161007_E_el1';


%% 1. Load data

cd 'E:\LaCie\DATA_KD\161007_E'

%load WorkspaceForProcessing_180705_evp001_test.mat
load BRFSWorkspace001_161007_E_brfs001ns2_2018-11-26.mat

Figuresdir = 'E:\LaCie\DATA_KD\BRFS_figs';
% 2. Plot Controll - Simultaneous Binocular Dioptic Grating 

%Create EV.A subsets - "biOriX" and "biOriY"
clear countX countY
a = 0;
countY = 0;
% % % % % for i = 1:length(Cond.A)
% % % % %     if Cond.A(i,3:8) == [60,60,1,1,0,2]
% % % % %         a = a+1;
% % % % %         biOriX(a,:) = EV.A(i,:);
% % % % %     elseif Cond.A(i,3:8) == [150,150,1,1,0,2]
% % % % %         countY = countY+1;
% % % % %         biOriY(countY,:) = EV.A(i,:);
% % % % %     end
% % % % % end

for i = 1:length(Cond.A)
    if Cond.A(i,3:8) == [75,75,.9,.9,0,2]
        a = a+1;
        biOriX(a,:) = EV.A(i,:);
    elseif Cond.A(i,3:8) == [165,165,.9,.9,0,2]
        countY = countY+1;
        biOriY(countY,:) = EV.A(i,:);
    end
end

% Calculate CSD for controll
clear stimLFP
pre = 100;
post = 500;

% Trigger LFP to stim on
for j = 1:length(biOriX)
    stimtmX = round(biOriX(j,1)/30); %divide by 30 to convert timepoints to 1kHz. LFP already in 1kHz
    refwinX = stimtmX-pre:stimtmX+post;
    stimLFPX(j,:,:) = LFP(refwinX,:); 
end

for k = 1:length(biOriY)
    stimtmY = round(biOriY(k,1)/30); %divide by 30 to convert timepoints to 1kHz. LFP already in 1kHz
    refwinY = stimtmY-pre:stimtmY+post;
    stimLFPY(k,:,:) = LFP(refwinY,:); 
end

% Compute average across trials
avgLFPX = squeeze(mean(stimLFPX,1));
avgLFPY = squeeze(mean(stimLFPY,1));

% Calculate CSD and subtract baseline
CSDX = calcCSD(avgLFPX(:,1:24)).*0.4;
CSDY = calcCSD(avgLFPY(:,1:24)).*0.4;
bl =1:50;
bl_CSDX = mean(CSDX(:,bl),2);
bl_CSDY = mean(CSDY(:,bl),2);
blCSDX = (CSDX-bl_CSDX); 
blCSDY = (CSDY-bl_CSDY);

% FIX

 %%%%  PLOT THE STANDARD ERROR ON THE MEAN ?

% 

% Plot controll
chan = Sink -1; %this ACTUALLY gets to the sink. Deepest sink found on this chan for 161005
tvec = -pre:post;
figure(1),cla
	h1= plot(tvec,blCSDX(chan,:)); hold on
	h2 = plot(tvec,blCSDY(chan,:));
        vline(0); 
        title('161007 E 001 75vs165 Pref ori CSD sink lines');        
        ylabel( 'nA/(mm^3)') ; 
        xlabel('time (ms)');
        legend([h1 h2],{'75','165'},'location','best');
        hold off
        
 saveas(gcf,fullfile(Figuresdir,strcat(Session,'_oriTune')), 'jpeg');

clearvars -except LFP EV Cond Sink Session Figuresdir
        
%% 3. Figure set 2, Xdeg on flash, (second presentation) CSD

pre = 100;
post = 500;

%Create EV.D subset - dCOSX
a = 0;
for i = 1:length(Cond.D)
    if Cond.D(i,3:8) == [165,75,.9,.9,800,3]
        a = a+1;
        dCOSoriX(a,:) = EV.D(i,:);
    end
end

% Trigger LFP to stim on
for j = 1:length(dCOSoriX)
    stimtm_dCOSoriX = round(dCOSoriX(j,1)/30); %divide by 30 to convert timepoints to 1kHz. LFP already in 1kHz
    refwin_dCOSoriX = stimtm_dCOSoriX-pre:stimtm_dCOSoriX+post;
    stimLFP_dCOSoriX(j,:,:) = LFP(refwin_dCOSoriX,:); 
end
 
% Compute average across trials
avgLFP_dCOSoriX = squeeze(mean(stimLFP_dCOSoriX,1));

% Calculate CSD and subtract baseline
CSD_dCOSoriX = calcCSD(avgLFP_dCOSoriX(:,1:24)).*0.4;
bl =1:50;
bl_CSD_dCOSoriX = mean(CSD_dCOSoriX(:,bl),2);
blCSD_dCOSoriX = (CSD_dCOSoriX-bl_CSD_dCOSoriX); 

% pad 
padCSD_dCOSoriX = padarray(blCSD_dCOSoriX,[1 0],NaN,'replicate'); %48x1301

% pull out line plots for later
    %supragranular 
    supGranChan = Sink - 5;
    supGran_X = padCSD_dCOSoriX(supGranChan,:);
    
    %granular
    granChan = Sink;
    Gran_X = padCSD_dCOSoriX(granChan,:);
    
    %infragranular
    infraGranChan = Sink + 5;
    infraGran_X = padCSD_dCOSoriX(infraGranChan,:);
    
% filter
gauss_sigma = 0.1;
fCSD_dCOSoriX = filterCSD(padCSD_dCOSoriX,gauss_sigma);

% plot filtered data
tvec = (-pre:post);
chanvec_dCOSoriX = linspace(1,24,size(fCSD_dCOSoriX,1));

figure(2)

    imagesc(tvec,chanvec_dCOSoriX,fCSD_dCOSoriX);
        colormap(flipud(jet)); colorbar;  vline(0);
        set(gca,'XTickLabel',[700 800 900 1000 1100 1200 1300]);
        title('161007 E 001 el 1 dCOS 75 presented on flash. non pref stim on flash');
        xlabel('time (ms)') 
        ylabel('electrode number') 
        mn1 = min(min(fCSD_dCOSoriX)); mx1 = max(max(fCSD_dCOSoriX));
        maxvallX = max([abs(mn1) abs(mx1)]);
        caxis([-maxvallX maxvallX]);
% %         caxis([-2390 1853]);

 saveas(gcf,fullfile(Figuresdir,strcat(Session,'_filtX')), 'jpeg');        
        
figure(3)        
    ax1 = subplot(3,1,1);
    	plot(supGran_X);
        ylim([-4500 1500])
        xlim([0 600])
        set(gca,'XTickLabel',[700 800 900 1000 1100 1200 1300]);
        title('dCOS 50 supragranular');
    ax2 = subplot(3,1,2);
        plot(Gran_X);
        ylim([-4500 1500]); xlim([0 600])
        set(gca,'XTickLabel',[700 800 900 1000 1100 1200 1300]);
        title('dCOS 50 Granular. Min = -1624.4');
    ax3 = subplot(3,1,3);
        plot(infraGran_X)
        ylim([-4500 1500]); xlim([0 600])
        set(gca,'XTickLabel',[700 800 900 1000 1100 1200 1300]);
        title('dCOS 50 infragranular');
        
 saveas(gcf,fullfile(Figuresdir,strcat(Session,'_lineX')), 'jpeg');

       clearvars -except LFP EV Cond Sink fCSD_dCOSoriX Session Figuresdir Gran_X maxvallX
        
%% 4. Figure set 3, Ydeg on flash, (second presentation) CSD

pre = 100;
post = 500;

%Create EV.D subset - dCOSX
a = 0;
for i = 1:length(Cond.D)
    if Cond.D(i,3:8) == [75,165,.9,.9,800,3]
        a = a+1;
        dCOSoriY(a,:) = EV.D(i,:);
    end
end

% Trigger LFP to stim on
for j = 1:length(dCOSoriY)
    stimtm_dCOSoriY = round(dCOSoriY(j,1)/30); %divide by 30 to convert timepoints to 1kHz. LFP already in 1kHz
    refwin_dCOSoriY = stimtm_dCOSoriY-pre:stimtm_dCOSoriY+post;
    stimLFP_dCOSoriY(j,:,:) = LFP(refwin_dCOSoriY,:); 
end
 
% Compute average across trials
avgLFP_dCOSoriY = squeeze(mean(stimLFP_dCOSoriY,1));

% Calculate CSD and subtract baseline
CSD_dCOSoriY = calcCSD(avgLFP_dCOSoriY(:,1:24)).*0.4;
bl =1:50;
bl_CSD_dCOSoriY = mean(CSD_dCOSoriY(:,bl),2);
blCSD_dCOSoriY = (CSD_dCOSoriY-bl_CSD_dCOSoriY); 

% pad 
padCSD_dCOSoriY = padarray(blCSD_dCOSoriY,[1 0],NaN,'replicate'); %48x1301

% pull out line plots for later
    %supragranular 
    supGranChan = Sink - 5;
    supGran_Y = padCSD_dCOSoriY(supGranChan,:);
    
    %granular
    granChan = Sink;
    Gran_Y = padCSD_dCOSoriY(granChan,:);
    
    %infragranular
    infraGranChan = Sink + 5;
    infraGran_Y = padCSD_dCOSoriY(infraGranChan,:);
    
% filter
gauss_sigma = 0.1;
fCSD_dCOSoriY = filterCSD(padCSD_dCOSoriY,gauss_sigma);

% plot filtered data
tvec = (-pre:post);
chanvec_dCOSoriY = linspace(1,24,size(fCSD_dCOSoriY,1));


figure(4)

    imagesc(tvec,chanvec_dCOSoriY,fCSD_dCOSoriY);
        colormap(flipud(jet)); colorbar;  vline(0);
        set(gca,'XTickLabel',[700 800 900 1000 1100 1200 1300]);
        title('161007 E 001 el 1 dCOS 165deg on flash.  pref stim on flash');
        xlabel('time (ms)') 
        ylabel('electrode number') 
        mn1 = min(min(fCSD_dCOSoriY)); mx1 = max(max(fCSD_dCOSoriY));
        maxvallY = max([abs(mn1) abs(mx1)]);
% % %          caxis([-maxvallY maxvallY]);
         caxis([-maxvallX maxvallX])

 saveas(gcf,fullfile(Figuresdir,strcat(Session,'_filtY')), 'jpeg');

figure(5)        
    ax1 = subplot(3,1,1);
    	plot(supGran_Y);
        ylim([-4500 1500]); xlim([0 600])
        set(gca,'XTickLabel',[700 800 900 1000 1100 1200 1300]);
        title('dCOS 165 supragranular');
    ax2 = subplot(3,1,2);
        plot(Gran_Y);
        ylim([-4500 1500]); xlim([0 600])
        set(gca,'XTickLabel',[700 800 900 1000 1100 1200 1300]);
        title('dCOS 165 Granular. Min = -796.2');
    ax3 = subplot(3,1,3);
        plot(infraGran_Y)
        ylim([-4500 1500]); xlim([0 600])
        set(gca,'XTickLabel',[700 800 900 1000 1100 1200 1300]);
        title('dCOS 165 infragranular');

saveas(gcf,fullfile(Figuresdir,strcat(Session,'_lineY')), 'jpeg');  

mnX = min(Gran_X);
mnY = min(Gran_Y);
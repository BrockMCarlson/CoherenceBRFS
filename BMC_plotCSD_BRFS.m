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

% COMPUTE AVERAGE ACROSS TRIALS
avgLFP = squeeze(mean(stimLFP,1));
avgLFP_elC = avgLFP(:,1:32);

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
fCSD_elC = filterCSD(padCSD_elC(2:23,:),gauss_sigma);
fCSD_elD = filterCSD(padCSD_elD(2:24,:),gauss_sigma);

%% Plot CSD
tvec = (-pre:post);
chanvec_elC = linspace(2,23,size(fCSD_elC,2));
chanvec_elD = linspace(2,24,size(fCSD_elD,2));

figure(3)
subplot(1,2,1); 
imagesc(tvec,chanvec_elC,fCSD_elC);
colormap(flipud(jet)); colorbar; vline(0);
mn1 = min(min(fCSD_elC)); mx1 = max(max(fCSD_elC));
maxval1 = max([abs(mn1) abs(mx1)]);
caxis([-maxval1 maxval1]);
title('161007 sinkbtm 15');

subplot(1,2,2);
imagesc(tvec,chanvec_elD,fCSD_elD);
colormap(flipud(jet)); colorbar; vline(0);
mn2 = min(min(fCSD_elD)); mx2 = max(max(fCSD_elD)); 
maxval2 = max([abs(mn2) abs(mx2)]);
caxis([-maxval2 maxval2]);
title('161007 sinkbtm 16');



dateFormatOut = 'yyyy-mm-dd';
saveDate = datestr(now,dateFormatOut);
saveName = strcat('fig_CSD_',Filename,'_',saveDate);
clear dateFormatOut saveDate
saveas(figure(3),saveName)


%CONGRATS... it finished
load gong
sound(y,Fs)



%% Helper functions
function CSD = calcCSD(avgLFP)
totchan = size(avgLFP,2);
if ndims(avgLFP) == 3
    maxtr = size(avgLFP,3);
else
    maxtr = 1;
end
%%% SET CSD PARAMS %%%
% electrode contact spacing in mm:
el_pos = [0.1:0.1:totchan/10];
N = length(el_pos); d = mean(diff(el_pos));
for i=1:N-2
    for j=1:N
        if (i == j-1)
            out(i,j) = -2/d^2;
        elseif (abs(i-j+1) == 1)
            out(i,j) = 1/d^2;
        else
            out(i,j) = 0;
        end
    end
end
%%%%%% END of CSD Params %%%%%%
for tr=1:maxtr
    CSD(:,:,tr)=(-1*out*squeeze(avgLFP(:,:,tr))');
end

end


%filterCSD.m
function [CSDf] = filterCSD(CSD, gauss_sigma)

if nargin < 2
    gauss_sigma = 0.1;
end

    new_CSD_matrix=[];
    totchan = (size(CSD,1) + 2)/10;
    el_pos = .1:.1:totchan;
    %npoints = 200; % number of points to plot in the vertical direction
    npoints = 10* size(CSD,1);
    le = length(el_pos)-2;
    first_z = el_pos(1)-(el_pos(2)-el_pos(1))/2; %plot starts at z1-h/2;
    last_z = el_pos(le)+(el_pos(le)-el_pos(le-1))/2; %ends at zN+h/2;
    zs = first_z:(last_z-first_z)/npoints:last_z;
    el_pos(le+1) = el_pos(le)+(el_pos(le)-el_pos(le-1)); % need an extra pos in for-loop
    j=1; %counter
    for i=1:length(zs) % all new positions
        if zs(i)>(el_pos(j)+(el_pos(j+1)-el_pos(j))/2) % > el_pos(j) + h/2
            j = min(j+1,le);
        end
        new_CSD_matrix(i,:) = CSD(j,:);
    end
    
    % and filter
    filter_range = 5*gauss_sigma; % numeric filter must be finite in extent
    [zs,CSDf]  = gaussian_filtering(zs,new_CSD_matrix,gauss_sigma,filter_range);
    
end


function [new_positions,gfiltered_CSD] = gaussian_filtering(positions,unfiltered_CSD,gauss_sigma,filter_range)
%function [new_positions, gfiltered_CSD]= ...
%gaussian_filtering(positions,unfiltered_CSD,gauss_sigma,filter_range)
%
%This function filters the CSD using a gaussian filter.
%
%positions: The CSD positions
%unfiltered_CSD: The unfiltered CSD matrix
%gauss_sigma: standard deviation of the gaussian
%filter_range: the filter width, default: 5*gauss_sigma

%Copyright 2005 Klas H. Pettersen under the General Public License,
%
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or any later version.
%
%See: http://www.gnu.org/copyleft/gpl.html

if nargin<4; filter_range = 5*gauss_sigma; end;
%first_position = positions(1);
%last_position = positions(length(positions));
%errormsg('CSD_matrix must have rows equal to length of positions.')

step = positions(2)-positions(1);
filter_positions = -filter_range/2:step:filter_range/2;
gaussian_filter = 1/(gauss_sigma*sqrt(2*pi))*exp(-filter_positions.^2/(2*gauss_sigma^2));
% figure();
% plot(filter_positions,gaussian_filter);
filter_length = length(gaussian_filter);
[m,n]=size(unfiltered_CSD);
temp_CSD=zeros(m+2*filter_length,n);
temp_CSD(filter_length+1:filter_length+m,:)=unfiltered_CSD(:,:); % one filter length of zeros on each side
scaling_factor = sum(gaussian_filter);
temp_CSD = filter(gaussian_filter/scaling_factor,1,temp_CSD); % filter works such that the first filter_length positions is crap
%gfiltered_CSD=temp_CSD(filter_length+1:2*filter_length-1+m,:);
gfiltered_CSD=temp_CSD(round(1.5*filter_length)+1:round(1.5*filter_length)+m,:); % first filter_length is crap, next 0.5 filter length corresponds to positions smaller than the original positions
%first_position = positions(1)-(filter_length-1)/2*step;
%last_position = positions(length(positions))+(filter_length-1)/2*step;
%new_positions = first_position:step:last_position;
new_positions = positions;
end
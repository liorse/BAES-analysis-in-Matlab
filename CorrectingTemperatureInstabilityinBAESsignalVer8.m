% Bug correction:
% 16 May 2013 - Reading Temperature directly from core data file. 
% 27 Dec 2012 - Taking real time measurements in between measurements and the
%               corresponding temperature value at these measurements time
% 31 Dec 2012 - Cross correlation is now fixing with respect to the
%               calibration ambient temperature
% 1 Jan 2014 - fix for reading when taking into account uncertainty data

% OUTPUT
% DeltaCavity how much did we have to move our wave in nm per measurment
% CavityCountsCorrected{1} for 365 corrected wavelength according to
% temperature changes
% CavityCountsCorrected{2} for 405

%% Read the BAES data files and Temperature log file
% filename input variable to the function
% input
FilenameData='20140810 stability test old filters.txt'; 
%FilenameTemp='Temp Exp2 day1.log';
NumberofHeaderLines=15; % may change in the future
deltaWavelength=0.01; % interpolate 0.01 nm
LinearFit = cell(2,1);
LinearFit{1} = [0.11709 -2.5525];  % Calibration for 365 
LinearFit{2} = [0.10515 -2.3030];  % Calibration for 405 

% get wavelength data - x coordinates from file
Data = importdata(FilenameData,',');
WaveLength = Data.data;

% Get counts data from 365 and 405 cavities
Data = importdata(FilenameData,',',NumberofHeaderLines);
Counts=Data.data(:,7:end);
CavityCounts=cell(2,1);
% replace this with a shift signal
ArtificialFlag = 0;
CavityCounts{1} = Counts(1:4:length(Counts(:,1)),:);
CavityCounts{2} = Counts(3:4:length(Counts(:,1)),:);
%CavityCounts{1} = CavityCountsArtificialData{1};
%CavityCounts{2} = CavityCountsArtificialData{2};


%% updated code 

% Interpolate temperature data to fit BAES data
EndTime = length(CavityCounts{1}(:,1));
%Timeline=ExperimentTime*(1:EndTime)./EndTime;
%CountsInterp=ExperimentTime*(1:length(Temp.data))./length(Temp.data);
%TempData=interp1(CountsInterp, Temp.data, Timeline,'spline');
Temperature=Data.data(1:4:length(Counts(:,1)),1);
OpticalPower=Data.data(1:4:length(Counts(:,1)),5);

Temperature(find(isnan(Temperature)==1))=mean(Temperature(1:2));
DeltaIndex = cell(2,1);
DeltaIndex{1} = round((LinearFit{1}(1)*Temperature+LinearFit{1}(2))./deltaWavelength);
DeltaIndex{2} = round((LinearFit{2}(1)*Temperature+LinearFit{2}(2))./deltaWavelength);

%% Fix the wavelength shift according to temperature variations
% 1. interpolate data for 0.01 nm resolution using splines between adjucent
% points
% 2. Use Temperature and calibration slope and offset calculated before to
% find appropriate deltas
% 3. Adjust the data wavelength accordingly
% 4. Interpolate back and draw corrected vs uncorrected and compare to
% interpolated data

% 1. interpolate data for 0.01 nm resolution using splines between adjucent
% points for all time points

freqRange=((WaveLength(end)-WaveLength(1))/round((WaveLength(end) - WaveLength(1))/deltaWavelength))*((1:(round((WaveLength(end) - WaveLength(1))/deltaWavelength)+1))-1) + WaveLength(1);
CavityCountsCorrected = cell(2,1); CavityCountsCorrected{1}= zeros(EndTime, length(WaveLength));CavityCountsCorrected{2}= zeros(EndTime, length(WaveLength));
CavityCountsCorrectedXcorr = cell(2,1); CavityCountsCorrectedXcorr{1} = zeros(EndTime, length(WaveLength));CavityCountsCorrectedXcorr{2} = zeros(EndTime, length(WaveLength));
CavityCountsCorrectedXcorrBGcorrected = cell(2,1); CavityCountsCorrectedXcorrBGcorrected{1} = zeros(EndTime, length(WaveLength));CavityCountsCorrectedXcorrBGcorrected{2} = zeros(EndTime, length(WaveLength));
CavityCountsWL_BG_Intensity = cell(2,1); CavityCountsWL_BG_Intensity{1} = zeros(EndTime, length(WaveLength));CavityCountsWL_BG_Intensity{2} = zeros(EndTime, length(WaveLength));
CavityCountsCorrectedBGcorrected = cell(2,1); CavityCountsCorrectedBGcorrected{1}= zeros(EndTime, length(WaveLength));CavityCountsCorrectedBGcorrected{2}=zeros(EndTime, length(WaveLength));
DeltaCavity = cell(2,1); DeltaCavity{1}=zeros(1,EndTime); DeltaCavity{2}=zeros(1,EndTime);
DeltaCavitynm = cell(2,1); DeltaCavitynm{1}=zeros(1,EndTime); DeltaCavitynm{2}=zeros(1,EndTime);

if ArtificialFlag == 1
    CavityCountsArtificialData = cell(2,1);
    CavityCountsArtificialData{1}= zeros(EndTime,length(WaveLength));CavityCountsArtificialData{2}= zeros(EndTime,length(WaveLength));
end

startIndex=1;



% correct for Intensity variations
%for j=startIndex:EndTime
%    d = sum(CavityCountsCorrectedXcorrBGcorrected{2}(j,:),2)./sum(CavityCountsCorrectedXcorrBGcorrected{2}(startIndex,:));
%    %d = CavityCountsCorrectedXcorrBGcorrected{2}(j,591)/CavityCountsCorrectedXcorrBGcorrected{2}(startIndex,591);
%    CavityCounts{2}(j,:) = CavityCountsCorrectedXcorrBGcorrected{2}(j,:)./d;
%end
DatesTime=datenum(Data.textdata(16:4:end,3));
ExperimentTime=[0; cumsum(DatesTime(2:end)-DatesTime(1:(length(DatesTime)-1)))*24];
for Cavity=1:2
    InterpResult=[];
    for SignalNum=startIndex:EndTime
        
        % interpolate data to 0.01 nm resolution for the specific cavity
        CountsInterp=interp1(WaveLength,CavityCounts{Cavity}(SignalNum,:), freqRange, 'spline');
        
        
        % fix the data according to temeprature and calibration
        if DeltaIndex{Cavity}(SignalNum)>0
            CountsInterpCorrected=[CountsInterp((abs(DeltaIndex{Cavity}(SignalNum))+1):end) zeros(abs(DeltaIndex{Cavity}(SignalNum)),1)'];
        elseif DeltaIndex{Cavity}(SignalNum)==0
            CountsInterpCorrected = CountsInterp;
        else
            CountsInterpCorrected=[zeros(abs(DeltaIndex{Cavity}(SignalNum)),1)' CountsInterp(1:end-(abs(DeltaIndex{Cavity}(SignalNum))))];
        end
        
        %interpolate back to 0.5 nm resolution
        CavityCountsCorrected{Cavity}(SignalNum,:)= interp1(freqRange, CountsInterpCorrected, WaveLength);
        
        % Fix the data according to deltas from cross correlation alone!
        
        if SignalNum==1
            % Save the the first signal as reference
            CountsInterpRef = CountsInterp;  % used to be CountsInterpCorrected
            DeltaCavitynm{Cavity}(SignalNum) = 0;
            DeltaCavity{Cavity}(SignalNum) =0;
        else
            % Prefrom a cross correlation calculation
             g=xcorr(CountsInterp,CountsInterpRef);
             DeltaCavitynm{Cavity}(SignalNum)= (find(g==max(g))-length(freqRange))*deltaWavelength;
             DeltaCavity{Cavity}(SignalNum) = find(g==max(g))-length(freqRange);
            % just a quick test
            
            %g=xcorr(CountsInterp(1:40000),CountsInterpRef(1:40000));
            %DeltaCavitynm{Cavity}(SignalNum)= (find(g==max(g))-length(freqRange(1:40000)))*deltaWavelength;
            %DeltaCavity{Cavity}(SignalNum) = find(g==max(g))-length(freqRange(1:40000));
            
        end
        
        
        if SignalNum==1
            % Save the the first signal as reference
            CavityCountsCorrectedXcorr{Cavity} = interp1(freqRange, CountsInterpCorrected, WaveLength); % used to be: CavityCounts{Cavity}(SignalNum,:)
        else
            sig=SignalNum; % 
            if DeltaCavity{Cavity}(sig)>0
                CountsInterpCorrected=[CountsInterp((abs(DeltaCavity{Cavity}(sig))+1):end) zeros(abs(DeltaCavity{Cavity}(sig)),1)'];
            elseif DeltaCavity{Cavity}(sig)==0
                CountsInterpCorrected = CountsInterp;
            else
                CountsInterpCorrected=[zeros(abs(DeltaCavity{Cavity}(sig)),1)' CountsInterp(1:end-(abs(DeltaCavity{Cavity}(sig))))];
            end
            
            %interpolate back to 0.5 nm resolution
            CavityCountsCorrectedXcorr{Cavity}(SignalNum,:) =interp1(freqRange, CountsInterpCorrected, WaveLength);
            
        end
        
        if ArtificialFlag == 1
            % Create artificial data according to the above temperature profile
            % and the specified slope. I should get the exact same slope back
            % and a perfect correction.
           
            if SignalNum==1
                % Save the the first signal as reference
                %CavityCountsArtificialData{Cavity} = CavityCounts{Cavity}(SignalNum,:); % ref
                CavityCounts{Cavity}(SignalNum,1:540)=zeros(1,540);
                CavityCounts{Cavity}(SignalNum,641:1024)=zeros(1,384);
                CountsInterp=interp1(WaveLength,CavityCounts{Cavity}(SignalNum,:), freqRange, 'pchip');
                
                CavityCountsArtificialDataInterp = CountsInterp;
            end    
            
            if DeltaIndex{Cavity}(SignalNum)>0
                CountsInterpCorrected=[CavityCountsArtificialDataInterp((abs(DeltaIndex{Cavity}(SignalNum))+1):end) zeros(abs(DeltaIndex{Cavity}(SignalNum)),1)'];
            elseif DeltaIndex{Cavity}(SignalNum)==0
                CountsInterpCorrected = CavityCountsArtificialDataInterp;
            else
                CountsInterpCorrected=[zeros(abs(DeltaIndex{Cavity}(SignalNum)),1)' CavityCountsArtificialDataInterp(1:end-(abs(DeltaIndex{Cavity}(SignalNum))))];
            end
            
            %interpolate back to 0.5 nm resolution
            CavityCountsArtificialData{Cavity}(SignalNum,:)= interp1(freqRange, CountsInterpCorrected, WaveLength, 'pchip');
            
            
            
        end
        
        if mod(SignalNum,100) == 0
            SignalNum
        end
        
    end
    
end

%% plot certain wave length as a function of time for the 405

FrequencyRange=395:0.3:426;
%FrequencyRange=405;
figure(1)

% plot uncorrected data
hold off

startIndex=1;

% Calculate a signals levels in reference with zero level noise at
% wavelengths far from 405 and 365
%for j=startIndex:EndTime
%    CavityCountsCorrectedBGcorrected{2}(j,:) = CavityCounts{2}(j,:) - (mean(CavityCounts{2}(j,800:900)));
%end

% Calculate a signals levels in reference with zero level noise at
% wavelengths far from 405 and 365
for j=startIndex:SignalNum
    CavityCountsCorrectedXcorrBGcorrected{2}(j,:) = CavityCountsCorrectedXcorr{2}(j,:) - (mean(CavityCountsCorrectedXcorr{2}(j,800:900)));
    CavityCountsCorrectedBGcorrected{2}(j,:) = CavityCountsCorrected{2}(j,:) - mean(CavityCountsCorrected{2}(j,800:900));
end

% correct for Intensity variations
for j=startIndex:SignalNum
    d = sum(CavityCountsCorrectedXcorrBGcorrected{2}(j,:),2)./sum(CavityCountsCorrectedXcorrBGcorrected{2}(startIndex,:));
    %d = CavityCountsCorrectedXcorrBGcorrected{2}(j,591)/CavityCountsCorrectedXcorrBGcorrected{2}(startIndex,591);
    CavityCountsWL_BG_Intensity{2}(j,:) = CavityCountsCorrectedXcorrBGcorrected{2}(j,:)./d;
end


for Frequency=FrequencyRange;
    
    FreqIndex=max(find(WaveLength(:)<=Frequency));
    
   
    plot(ExperimentTime,100*(CavityCounts{2}(startIndex:end,FreqIndex)-CavityCounts{2}(startIndex,FreqIndex))./CavityCounts{2}(startIndex,FreqIndex))
    hold all
   
end


legend('405','Location','Northwest')
ylabel('Change in signal(%)')
xlabel('Time (h)')
a=axis;
miny=a(3);
maxy=a(4)
axis([0 ExperimentTime(end) miny maxy])



%% plot certain wave length as a function of time for the 365

%FrequencyRange=360:0.3:395;
FrequencyRange=360:0.3:380;
figure(2)

% plot uncorrected data

hold off


startIndex=1;

% Calculate a signals levels in reference with zero level noise at
% wavelengths far from 405 and 365
for j=startIndex:SignalNum
    CavityCountsCorrectedXcorrBGcorrected{1}(j,:) = CavityCountsCorrectedXcorr{1}(j,:) - mean(CavityCountsCorrectedXcorr{1}(j,800:900));
    CavityCountsCorrectedBGcorrected{1}(j,:) = CavityCountsCorrected{1}(j,:) - mean(CavityCountsCorrected{1}(j,800:900));
end

% correct for Intensity variations
%
for j=startIndex:SignalNum
    d = sum(CavityCountsCorrectedXcorrBGcorrected{1}(j,:),2)./sum(CavityCountsCorrectedXcorrBGcorrected{1}(startIndex,:));
   % d = CavityCountsCorrectedXcorrBGcorrected{1}(j,512)/CavityCountsCorrectedXcorrBGcorrected{1}(startIndex,512);
    CavityCountsWL_BG_Intensity{1}(j,:) = CavityCountsCorrectedXcorrBGcorrected{1}(j,:)./d;
end

for Frequency=FrequencyRange;
    
    FreqIndex=max(find(WaveLength(:)<=Frequency));
    plot(ExperimentTime,100*(CavityCounts{1}(startIndex:end,FreqIndex)-CavityCounts{1}(startIndex,FreqIndex))./CavityCounts{1}(startIndex,FreqIndex))
    hold all
    
    
end


legend('365','Location','Northeast')
ylabel('Change in signal(%)')
xlabel('Time (h)')
a=axis;
miny=a(3);
maxy=a(4)
axis([0 ExperimentTime(end) miny maxy])



%%
figure(3)
subplot(4,1,1)
startIndex=1;
% plot the Temperature variation
plot(ExperimentTime, Temperature)
%legend('Temperature as a function time','Location','North')
ylabel('Temperature(deg)')

grid on

subplot(4,1,2)
hold off

%FreqIndex=find(CavityCounts{1}(1,:)==max(CavityCounts{1}(1,:),[],2));
FreqIndex=max(find(WaveLength(:)<=365)); %446.7 for HeNe laser
plot(ExperimentTime, (CavityCounts{1}(:,FreqIndex)-CavityCounts{1}(1,FreqIndex))/CavityCounts{1}(1,FreqIndex)*100)  % 365

hold all
%FreqIndex=find(CavityCounts{2}(1,:)==max(CavityCounts{2}(1,:),[],2));
FreqIndex=max(find(WaveLength(:)<=405))
plot(ExperimentTime, (CavityCounts{2}(:,FreqIndex)-CavityCounts{2}(1,FreqIndex))/CavityCounts{2}(1,FreqIndex)*100)  % 365
%plot(ExperimentTime(startIndex:SignalNum), (OpticalPower(:)-OpticalPower(1))./OpticalPower(1)*100)

legend('365', '405','Location','Northeast')
ylabel('Counts change (%)')

grid on

subplot(4,1,3)
hold off
plot(ExperimentTime, DeltaCavitynm{1}(:))  % 365
hold all
plot(ExperimentTime, DeltaCavitynm{2}(:))  % 405
legend('365', '405','Location','Northeast')
ylabel('Shift in wavelength (nm)')

grid on

subplot(4,1,4)
hold off
plot(ExperimentTime(startIndex:SignalNum), (max(CavityCountsCorrectedXcorr{1}(:,:),[],2)-max(CavityCountsCorrectedXcorr{1}(1,:)))/max(CavityCountsCorrectedXcorr{1}(1,:))*100)  % 365
hold all
% integration over all wavelength
%plot(ExperimentTime(startIndex:SignalNum), sum(CavityCountsCorrectedXcorr{1}(startIndex:SignalNum,:),2)./sum(CavityCountsCorrectedXcorr{1}(startIndex,:)))  % 365
plot(ExperimentTime(startIndex:SignalNum), (max(CavityCountsCorrectedXcorr{2},[],2)-max(CavityCountsCorrectedXcorr{2}(1,:)))/max(CavityCountsCorrectedXcorr{2}(1,:))*100)  % 365
legend('365','405')
ylabel('Counts at maximum(change in %)')
xlabel('Time (h)')
grid on



%subplot(4,1,4)

%plot(ExperimentTime(startIndex:SignalNum),max(CavityCountsCorrectedXcorr{2},[],2))  % 405
% integration over all wavelength
%plot(ExperimentTime(startIndex:SignalNum), sum(CavityCountsCorrectedXcorr{2}(startIndex:SignalNum,:),2)./sum(CavityCountsCorrectedXcorr{2}(startIndex,:)))  % 405

%legend('365 Intensity change', '365 integration Intensity change','405 Intensity change', '405 integration Intensity change','Location','SouthWest')


% for j=startIndex:SignalNum
%     d = CavityCountsCorrectedXcorr{2}(j,591)/CavityCountsCorrectedXcorr{2}(startIndex,591);
%     CavityCountsCorrectedXcorr{2}(j,:) = CavityCountsCorrectedXcorr{2}(j,:)./d;
% end


%% Calibration output
figure(4)

startIndex=1;
EndIndex=EndTime;
% plot 365 vs Temperature and perform a linear fit
subplot(2,1,1)
hold off
plot(Temperature, DeltaCavitynm{1},'*')  % 365
hold all
plot(Temperature(startIndex:EndIndex), DeltaCavitynm{1}(startIndex:EndIndex),'r*')  % 365

legend('365 xcorr vs Temp','Location','NorthWest')
ylabel('Shift from Xcorr (nm)')
xlabel('Temperature(deg)')
grid on


% Perform a linear fit
%[cfun,gof,output] = fit(TempData(2:end)',DeltaCavity{1}(:), 'poly1') %
[cfun,gof,output] = fit(Temperature(startIndex:EndIndex),DeltaCavitynm{1}(startIndex:EndIndex)', 'poly1')
plot(Temperature(startIndex:EndIndex), feval(cfun,Temperature(startIndex:EndIndex)))

str = sprintf('y = %0.5f*X %0.5f \nR square = %0.5f',cfun.p1,cfun.p2,gof.rsquare)
text(Temperature(round(median(1:EndTime))),DeltaCavitynm{1}(round(median(1:EndTime)))-0.05, str,'HorizontalAlignment','left')

% plot 405 vs Temeprature
subplot(2,1,2)
hold off
plot(Temperature, DeltaCavitynm{2} ,'*')  % 405
hold all
startIndex=1;
plot(Temperature(startIndex:EndIndex), DeltaCavitynm{2}(startIndex:EndIndex),'r*')  % 365
legend('405 xcorr vs Temp','Location','NorthWest')
ylabel('Shift from Xcorr (nm)')
xlabel('Temperature(deg)')
grid on

% Perform a linear fit
[cfun,gof,output] = fit(Temperature(startIndex:EndIndex),DeltaCavitynm{2}(startIndex:EndIndex)', 'poly1')
hold all
%plot(TempData(2:end), feval(cfun,TempData(2:end)))
plot(Temperature(startIndex:EndIndex), feval(cfun,Temperature(startIndex:EndIndex)))
str = sprintf('y = %0.5f*X %0.5f \nR square = %0.5f',cfun.p1,cfun.p2,gof.rsquare)

text(Temperature(round(median(1:EndTime))),DeltaCavitynm{2}(round(median(1:EndTime)))-0.05, str,'HorizontalAlignment','left')

%% plot Intensity as a function of Temperature
figure(5)
startIndex=1; 
EndIndex=EndTime;

subplot(2,1,1)
hold off
IntensityShift365 = (max(CavityCountsCorrectedXcorr{1},[],2)-max(CavityCountsCorrectedXcorr{1}(1,:)))/max(CavityCountsCorrectedXcorr{1}(1,:))*100;
plot(Temperature(startIndex:EndIndex), IntensityShift365(startIndex:EndIndex) ,'*')  % 365
hold all

subplot(2,1,2)
hold off
IntensityShift405 = (max(CavityCountsCorrectedXcorr{2},[],2)-max(CavityCountsCorrectedXcorr{2}(1,:)))/max(CavityCountsCorrectedXcorr{2}(1,:))*100;
plot(Temperature(startIndex:EndIndex), IntensityShift405(startIndex:EndIndex),'*')  % 405
hold all



% Perform a linear fit
startIndex=7; 
EndIndex=EndTime-7; %length(TempData);

subplot(2,1,1)
plot(Temperature(startIndex:EndIndex), IntensityShift365(startIndex:EndIndex),'r*')  % 365

subplot(2,1,2)
plot(Temperature(startIndex:EndIndex), IntensityShift405(startIndex:EndIndex),'r*')  % 365

%[cfun,gof,output] = fit(TempData(2:end)',DeltaCavity{1}(:), 'poly1') %
subplot(2,1,1)

[cfun,gof,output] = fit(Temperature(startIndex:EndIndex), IntensityShift365(startIndex:EndIndex), 'poly1')
plot(Temperature(startIndex:EndIndex), feval(cfun,Temperature(startIndex:EndIndex)))

str = sprintf('y = %0.5f*X %0.5f \nR square = %0.5f',cfun.p1,cfun.p2,gof.rsquare)
text(Temperature(round(median(1:EndTime))),IntensityShift365(round(median(1:EndTime)))-0.01, str,'HorizontalAlignment','left')

legend('365','Location','SouthWest')
ylabel('Relative shift in intensity')
xlabel('Temperature(deg)')
grid on

% plot 405 vs Temeprature
subplot(2,1,2)
legend('405','Location','SouthWest')
ylabel('Relative shift in intensity')
xlabel('Temperature(deg)')
grid on

% Perform a linear fit
[cfun,gof,output] = fit(Temperature(startIndex:EndIndex),IntensityShift405(startIndex:EndIndex), 'poly1')
hold all
%plot(TempData(2:end), feval(cfun,TempData(2:end)))
plot(Temperature(startIndex:EndIndex), feval(cfun,Temperature(startIndex:EndIndex)))
str = sprintf('y = %0.5f*X %0.5f \nR square = %0.5f',cfun.p1,cfun.p2,gof.rsquare)

text(Temperature(round(median(1:EndTime))),IntensityShift405(round(median(1:EndTime)))-0.01, str,'HorizontalAlignment','left')



%% print to file
%print -f1 'figure1' -dpng
%print -f2 'figure2' -dpng
print -f3 'figure3' -dpdf
%print -f4 'figure4' -dpng
%print -f5 'figure5' -dpng

figure(3)
% calculating total optical power
(sum(CavityCounts{1}(107,:))-sum(CavityCounts{1}(1,:)))/sum(CavityCounts{1}(1,:))*100

%% optical data
figure(7)
subplot(3,1,1)
startIndex=1;
% plot the Temperature variation
plot(ExperimentTime, Temperature)
%legend('Temperature as a function time','Location','North')
xlabel('Experiment time(h)')
ylabel('Temperature(deg)')

grid on

subplot(3,1,2)
hold off

plot(ExperimentTime(startIndex:SignalNum), (OpticalPower(:)-OpticalPower(1))./OpticalPower(1)*100)
grid on
legend('Optical Power','Location','Northeast')
ylabel('Counts change (%)')
xlabel('Experiment time(h)')

subplot(3,1,3)
plot(Temperature,(OpticalPower(:)-OpticalPower(1))./OpticalPower(1)*100,'*')
grid on
xlabel('Temperature (deg)')
ylabel('Counts change (%)')
%%
% createfigureOpticalPower(ExperimentTime,Temperature,(OpticalPower(:)-OpticalPower(1))./OpticalPower(1)*100)
print  'figure1' -dpdf
print  'figure2' -dpdf
print  'figure3' -dpdf


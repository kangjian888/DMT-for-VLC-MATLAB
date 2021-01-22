% Starting data acquisition and BER calculation
%% Global Parameter
N= 128;
M = 6; %16QAM modulation scheme
Ncp = 32;
RepeatTime = 2;
FrameNum = 100;
SubCarrierNum = 59; % every OFDM frame has 52 subcarriers
LowPaddingNum = 0; % padding number of low frequency
HighPaddingNum = 0; % padding number of high frequency
Pilot_pos = [7,21,43,57]; % this is the position of pilot
ThresholdValue = 50;

%% Synchronization
i = 1; % index of input data
% Some registers used in Frame Syn
L1 = 32;
L2 = 70;
L3 = N;
L4 = N*(RepeatTime-1) + Ncp;
SynchronizationOffset = 5;% this parameter is the offset position of the beginning of the sequence. 
                                  % if this value is 0, Then the sequence beginning posision is the first sample after the repeaded time peak
% Parameter Check
if SynchronizationOffset > Ncp
    error('The SynchronizationOffset is invalid, Please adjust the parameters.');
end

% Getting Data Processing
for j = 1: 10
    try
        % -----------------------------------------------------------
        % SyncPoint 'SettingsApplied' - all the settings were applied
        % -----------------------------------------------------------
        % Arming the SCOPE for single acquisition
        % -----------------------------------------------------------
        scope.SetTimeoutMilliseconds(2000); % Acquisition timeout - set it higher than the acquisition time
        scope.Write('SING');
        % -----------------------------------------------------------
        % DUT_Generate_Signal() - in our case we use Probe compensation signal
        % where the trigger event (positive edge) is reoccuring
        % -----------------------------------------------------------
        fprintf('Waiting for the acquisition to finish... ');
        tic
        scope.QueryString('*OPC?'); % Using *OPC? query waits until the instrument finished the Acquisition
        toc
        scope.ErrorChecking(); % Error Checking after the acquisition is finished
        % -----------------------------------------------------------
        % SyncPoint 'AcquisitionFinished' - the results are ready
        % -----------------------------------------------------------
        % Fetching the waveform in ASCII format
        % -----------------------------------------------------------
        samplesCount = scope.QueryInteger('ACQ:POIN?'); % Query the expected samples count
        fprintf('Fetching waveform in ASCII format... ');
        tic
        waveformASC_POS = scope.QueryASCII_ListOfDoubles('FORM ASC;:CHAN1:DATA?', samplesCount); % samplesCount is the maximum allowed samples to read
        %waveformASC_NEG = scope.QueryASCII_ListOfDoubles('FORM ASC;:CHAN4:DATA?', samplesCount); % samplesCount is the maximum allowed samples to read
        toc
        fprintf('Samples count: %d\n', size(waveformASC_POS, 2));
        %fprintf('Samples count: %d\n', size(waveformASC_NEG, 2));
        scope.ErrorChecking(); % Error Checking after the data transfer
        % -----------------------------------------------------------
        % Fetching the trace in Binary format
        % Transfer of traces in binary format is faster.
        % The waveformBIN data and waveformASC data are however the same.
        % -----------------------------------------------------------
        %fprintf('Fetching waveform in binary format... ');
        %tic
        %waveformBIN = scope.QueryBinaryFloatData('FORM REAL,32;:CHAN1:DATA?');
        %toc
        %fprintf('Samples count: %d\n', size(waveformBIN, 2));
        %scope.ErrorChecking(); % Error Checking after the data transfer
        %plot(waveformBIN); % Displaying the waveform
        % -----------------------------------------------------------
        % Making an instrument screenshot and transferring the file to the PC
        % -----------------------------------------------------------
         scope.Write('CHAN1:STAT ON'); % Switch Channel 1 ON
         %scope.Write('CHAN4:STAT ON'); % Switch Channel 1 ON
        % -----------------------------------------------------------
        % Closing the session
        % -----------------------------------------------------------
        %scope.Close() % Closing the session to the instrument
        % -----------------------------------------------------------
        % Error handling
        % -----------------------------------------------------------
        waveformASC = waveformASC_POS;
        LengthReceivedSeq = length(waveformASC);
        DataIn = 0;
        FirstOrderShiftRam = zeros(1,L1);
        SecondOrderShiftRam = zeros(1,L2);
        CorrelatedSTSRam = zeros(1,L1);
        CorrelatedOriRam = zeros(1,L1);
        CorrelatedSum = 0;
        CorrelatedOriSum = 0;
        
        % Some registers used in Symbol Syn
        DataShiftRam = zeros(1,L3);
        SignDataShiftRam = zeros(1,L3);
        LongTSTRam = zeros(1,L3);
        SignLongTSTRam = zeros(1,L3);
        DataShiftForChannelEstimationRam = zeros(1,L4);
        SymbolDetected = 0;
        TimeCounter = 0;%detect the range between the peak.
        BitSumPro = 0;
        Synchronized = 0;
        
        % Some register used in removing CP module

        CpDataIn = 0;
        FrameDetcted = 0; % every time ratio is larger than threshold, this value will be added to 1.
        SymbolCounter = 0; % counting the symbol after the synchronization
        FrameCounter = 0;% counting the OFDM frame after the synchronization
        SymbolSynEnable = 0;
        
        % Some registers used in Channel estimation
        OneLSTafterFFT = zeros(1,N);
        AccumulateLTSafterFFT = zeros(1,N);
        
        % Some registers used in demodulation
        OFDMSymbolWithoutPadding = [];
        plot_enable = 0;
        i = 0;
        RatioOutputArray = [];
        MproArray = [];
        BitSumProAarry = [];
        BitSumAarry = [];


        while(1)
            i = i + 1;
            % the interpolation will be added later
            DataInMSB = DataIn;
            if i <= LengthReceivedSeq
                DataIn = waveformASC(i);
            else
                DataIn = 0; %When the data transmission is done, the 0 or noise is transmitted in the system
                i = i-1;
                break;
                disp("Symbol synchronization cannot be found in the whole seq ...");
            end
            FirstOrderShiftRamMSB = FirstOrderShiftRam(L1);
            FirstOrderShiftRam = [DataInMSB,FirstOrderShiftRam(1:L1-1)];
            SecondOrderShiftRamMSB = SecondOrderShiftRam(L2);
            SecondOrderShiftRam = [FirstOrderShiftRamMSB,SecondOrderShiftRam(1:L2-1)];
            CorrelatedSTS = sign(DataInMSB) * sign(FirstOrderShiftRamMSB);
            CorrelatedSTSRamMSB = CorrelatedSTSRam(L1);
            CorrelatedSTSRam = [CorrelatedSTS,CorrelatedSTSRam(1:L1-1)];
            CorrelatedSum = CorrelatedSum + CorrelatedSTS - CorrelatedSTSRamMSB;
            CorrelatedOri = sign(DataInMSB)*sign(DataInMSB);
            CorrelatedOriRamMSB = CorrelatedOriRam(L1);
            CorrelatedOriRam = [CorrelatedOri,CorrelatedOriRam(1:L1-1)];
            CorrelatedOriSum = CorrelatedOriSum + CorrelatedOri - CorrelatedOriRamMSB;
            % This part is used to plot the frame synchronization result
            if CorrelatedOriSum == 0
                RatioOutput = 0;
            else
                RatioOutput = abs(CorrelatedSum)/CorrelatedOriSum;
            end
            RatioOutputArray(i) = RatioOutput;
            % Below is Frame detection module
            if CorrelatedSum > CorrelatedOriSum/2 %threshold is 0.5
                FrameDetcted = FrameDetcted  + 1;
            end
            if FrameDetcted == 32
                disp('Frame is detected ...');
                SymbolSynEnable = 1;
            end
            % Symbol Synchronization
            if(SymbolSynEnable)
                DataShiftRamMSB = DataShiftRam(L3);
                DataShiftRam = [SecondOrderShiftRamMSB,DataShiftRam(1:L3-1)];
                DataShiftForChannelEstimationRamMSB = DataShiftForChannelEstimationRam(L4);
                DataShiftForChannelEstimationRam = [DataShiftRamMSB,DataShiftForChannelEstimationRam(1:L4-1)];
                [~,LongTSTRam,~]= TrainingSeqGenParkMethod(N,Ncp,1);
                LongTSTRam = fliplr(LongTSTRam);
                SignDataShiftRam(find(DataShiftRam>=0)) = 0;
                SignDataShiftRam(find(DataShiftRam<0)) = 1;
                SignLongTSTRam(find(LongTSTRam>=0)) = 0;
                SignLongTSTRam(find(LongTSTRam<0)) = 1;
                BitSum = sum(~xor(SignDataShiftRam,SignLongTSTRam));
                BitSumAarry(i) = BitSum;
                BitSumLastPro = BitSumPro;
                BitSumPro = 2*BitSum-(N);
                BitSumProAarry(i) = BitSumPro;
                if BitSumPro > 60 % first threshold
                    Mpro = BitSumPro - BitSumLastPro;
                else
                    Mpro = BitSumPro;
                end
                MproArray(i) = Mpro;
                if Mpro > ThresholdValue %Threshold is 120
                    if TimeCounter == 0 
                        SymbolDetected = SymbolDetected + 1;
                    elseif TimeCounter > N/2 % time range is larger than 
                        SymbolDetected = SymbolDetected + 1;
                        TimeCounter = 0;
                    end
                end
                if SymbolDetected ~= 0 && SymbolDetected ~= RepeatTime
                    TimeCounter = TimeCounter + 1;
                end
                if SymbolDetected == RepeatTime
                    Synchronized = 1;
                    SymbolDetected = 0;
                    disp('The Symbol is synchronized ...');
                    PayLoadLength = LengthReceivedSeq + L1 + L2 - i;
                end
                if Synchronized
                    if(FrameCounter <= RepeatTime-1)
                        SymbolCounter = SymbolCounter + 1;
                        CpDataIn(SymbolCounter) = DataShiftForChannelEstimationRam((RepeatTime-1)*N+SynchronizationOffset);
                        if SymbolCounter == N
                            FrameCounter = FrameCounter + 1;
                            SymbolCounter = 0;
                            OneLTSafterFFT = fft(CpDataIn); %FFT transformation
                            AccumulateLTSafterFFT = AccumulateLTSafterFFT + OneLTSafterFFT;
                            if FrameCounter == RepeatTime
                                AverageLTSafterFFT = AccumulateLTSafterFFT./RepeatTime;
                                ReceivedLTSFreqDom = AverageLTSafterFFT(2:N/2);
                                [~,~,LTSFreqDom] = TrainingSeqGenParkMethod(N,Ncp,RepeatTime);
                                H = ReceivedLTSFreqDom./LTSFreqDom(2:N/2); % Channel information 
                            end
                        end
                    else
                        if PayLoadLength-(FrameCounter-RepeatTime+1)*(N+Ncp) >= 0 && FrameCounter - RepeatTime < FrameNum
                            SymbolCounter = SymbolCounter + 1;
                            CpDataIn(SymbolCounter) = DataShiftForChannelEstimationRam((RepeatTime-1)*N+SynchronizationOffset);
                        else
                            disp('Received Done, Demodulation processing begin ...');
                            plot_enable = 1;
                            if FrameCounter-RepeatTime == 0
                                disp("the sequence is too short to include any effective payload.")
                                break;
                            end
                            EffectiveDataSymbolReshape = reshape(OFDMSymbolWithoutPadding,1,SubCarrierNum*(FrameCounter-RepeatTime));
                            % demodulation
                            DataSerial= qamdemod(EffectiveDataSymbolReshape.',2^M,'gray','OutputType','bit','PlotConstellation',false,'UnitAveragePower',true);
                             [~, DataSerialOri] = prbs_gen(7, 1, 1, length(DataSerial));
                            [bit_error_num, bit_error_rate] = biterr(DataSerial',DataSerialOri);
                            str = ["bit error rate is" ,num2str(bit_error_rate)];
                            disp(str);
                            % constellation calculation and MSE calculation
                            SymbolOnSpecificSubcarrier = OFDMSymbolWithoutPadding(4,:);
                            SymbolOnSpecificSubcarrier_2 = OFDMSymbolWithoutPadding(59,:);
                            SymbolOnallSubcarrier = reshape(OFDMSymbolWithoutPadding(1:20,:),1,20*size(OFDMSymbolWithoutPadding,2));
                            break;
                        end
                        if SymbolCounter == N+Ncp
                            FrameCounter = FrameCounter + 1;
                            SymbolCounter = 0;
                            % removing CP
                            OneOFDMSymbolWithoutCP = CpDataIn(Ncp+1:end);
                            % FFT transformation
                            OneOFDMSymbolFreqDom = fft(OneOFDMSymbolWithoutCP);
                            OneOFDMSymbolWithoutHermitian = OneOFDMSymbolFreqDom(2:N/2);
                            OneOFDMSymbolAfterChannelEstimation = OneOFDMSymbolWithoutHermitian./H;
                            % Pilot extraction
                            RecievedPilot = OneOFDMSymbolAfterChannelEstimation(Pilot_pos);
                            OneOFDMSymbolWithoutPilot = OneOFDMSymbolAfterChannelEstimation([1:Pilot_pos(1)-1,Pilot_pos(1)+1:Pilot_pos(2)-1,Pilot_pos(2)+1:Pilot_pos(3)-1,Pilot_pos(3)+1:Pilot_pos(4)-1,Pilot_pos(4)+1:N/2-1]);
                            OneOFDMSymbolWithoutPadding = OneOFDMSymbolWithoutPilot(LowPaddingNum+1:LowPaddingNum+SubCarrierNum);
                            OFDMSymbolWithoutPadding = [OFDMSymbolWithoutPadding,OneOFDMSymbolWithoutPadding.'];
                        end
                    end
                end 
            end        
        end
        figure(4)
        plot(1:i,RatioOutputArray);
        disp('Drawing the Symbol Synchronization ...');
        figure(5)
        plot(1:i,MproArray,'r');
        if (plot_enable)
            %plot(1:i,MproArray,'r',1:i,BitSumProAarry,'b',1:i,BitSumAarry,'g');
            figure(6)
            scatter(real(SymbolOnSpecificSubcarrier),imag(SymbolOnSpecificSubcarrier),'*');
            figure(7)
            scatter(real(SymbolOnSpecificSubcarrier_2),imag(SymbolOnSpecificSubcarrier_2),'*');
            figure(8)
            scatter(real(SymbolOnallSubcarrier),imag(SymbolOnallSubcarrier),'*');            
            plot_enable = 0;
        end
    catch ME
        switch ME.identifier
            case 'VISA_Instrument:ErrorChecking'
                % Perform your own additional steps here
                rethrow(ME);
            otherwise
                rethrow(ME);
        end
    end

end
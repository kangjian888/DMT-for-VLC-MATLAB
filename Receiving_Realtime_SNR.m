% Starting data acquisition and BER calculation
%% Global Parameter
N= 128;
M = 3; %16QAM modulation scheme
Ncp = 32;
RepeatTime = 2;
FrameNum = 200; % ofdm symbol number of each package transmitted
TotalFrameNum = 10000; % total ofdm symbol number you want to get
SubCarrierNum = 59; % every OFDM frame has 52 subcarriers
LowPaddingNum = 0; % padding number of low frequency
HighPaddingNum = 0; % padding number of high frequency
Pilot_pos = [7,21,43,57]; % this is the position of pilot
ThresholdValue = 60;

POLY_LENGTH = 7;
POLY_TAP = 1;
PILOT_NUM = 4; % there is four subcarriers are used for pilot
NBITS = 1; %Data width of parallel output
NUM = FrameNum*SubCarrierNum*M; %Byte number

WORD_LENGTH = 16;
FRACTION_LENGTH = 14;

%% Synchronization
i = 1; % index of input data
% Some registers used in Frame Syn
L1 = 32;
L2 = 70;
L3 = N;
L4 = N*(RepeatTime-1) + Ncp;
SynchronizationOffset = 6;% this parameter is the offset position of the beginning of the sequence. 
                                  % if this value is 0, Then the sequence beginning posision is the first sample after the repeaded time peak
SpecificSubcarrier = 40;
% Parameter Check
if SynchronizationOffset > Ncp
    error('The SynchronizationOffset is invalid, Please adjust the parameters.');
end

% globle registers and varibles
OFDMSymbolWithoutPaddingAcc = [];
OFDMSymbolNumArray = [];
% Getting Data Processing

[~,DataSymbolQuanReshape,DataSerialOri]= TransmitterTopNew(N,M,Ncp,RepeatTime,FrameNum,LowPaddingNum,HighPaddingNum,SubCarrierNum,Pilot_pos,WORD_LENGTH,FRACTION_LENGTH);
while(1)
    try
        scope.SetTimeoutMilliseconds(2000); % Acquisition timeout - set it higher than the acquisition time
        scope.Write('SING');
        fprintf('Waiting for the acquisition to finish... ');
        tic
        scope.QueryString('*OPC?'); % Using *OPC? query waits until the instrument finished the Acquisition
        toc
        scope.ErrorChecking(); % Error Checking after the acquisition is finished
        samplesCount = scope.QueryInteger('ACQ:POIN?'); % Query the expected samples count
        fprintf('Fetching waveform in ASCII format... ');
        tic
        waveformASC_TRAN = scope.QueryASCII_ListOfDoubles('FORM ASC;:CHAN3:DATA?', samplesCount); % samplesCount is the maximum allowed samples to read
        waveformASC_POS = scope.QueryASCII_ListOfDoubles('FORM ASC;:CHAN2:DATA?', samplesCount); % samplesCount is the maximum allowed samples to read
        waveformASC_NEG  = scope.QueryASCII_ListOfDoubles('FORM ASC;:CHAN1:DATA?', samplesCount); % samplesCount is the maximum allowed samples to read
        toc
        fprintf('Samples count: %d\n', size(waveformASC_POS, 2));
        fprintf('Samples count: %d\n', size(waveformASC_NEG, 2));
        scope.ErrorChecking(); % Error Checking after the data transfer
        scope.Write('CHAN1:STAT ON'); % Switch Channel 1 ON
        scope.Write('CHAN3:STAT ON'); % Switch Channel 1 ON
        scope.Write('CHAN2:STAT ON'); % Switch Channel 1 ON
        %waveformASC = waveformASC_TRAN;
        %waveformASC = waveformASC_POS;
        waveformASC = waveformASC_POS-waveformASC_NEG;
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
        i = 1;
        RatioOutputArray = [];
        MproArray = [];
        BitSumProAarry = [];
        BitSumAarry = [];


        while(1)
            % the interpolation will be added later
            DataInMSB = DataIn;
            if i <= LengthReceivedSeq
                DataIn = waveformASC(i);
            elseif i > LengthReceivedSeq && i <= LengthReceivedSeq + 1 + L1 + L2 + L3 +L4
                DataIn = 0; %When the data transmission is done, the 0 or noise is transmitted in the system
                %i = i-1;
            else
                disp("Symbol synchronization cannot be found in the whole seq ...");
                break;
            end
            i = i + 1;
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
                    PayLoadLength = LengthReceivedSeq + L1 + L2 - i + 1;
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
                            OFDMSymbolWithoutPaddingAcc = [OFDMSymbolWithoutPaddingAcc,OFDMSymbolWithoutPadding]; %OFDM symbols smapled in several times
                            % demodulation
                            OFDMSymbolNumArray = [OFDMSymbolNumArray,FrameCounter-RepeatTime]; % the OFDM number getted each sampling time
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
                            % Pilot extraction, no channel estimation
                            OneOFDMSymbolWithoutPilot = OneOFDMSymbolAfterChannelEstimation([1:Pilot_pos(1)-1,Pilot_pos(1)+1:Pilot_pos(2)-1,Pilot_pos(2)+1:Pilot_pos(3)-1,Pilot_pos(3)+1:Pilot_pos(4)-1,Pilot_pos(4)+1:N/2-1]);
                            OneOFDMSymbolWithoutPadding = OneOFDMSymbolWithoutPilot(LowPaddingNum+1:LowPaddingNum+SubCarrierNum);
                            OFDMSymbolWithoutPadding = [OFDMSymbolWithoutPadding,OneOFDMSymbolWithoutPadding.'];
                        end
                    end
                end 
            end        
        end
        disp('Drawing the Symbol Synchronization ...');
        figure(1)
        plot(1:i,RatioOutputArray);
        figure(2)
        plot(1:i,MproArray,'r');
    catch ME
        switch ME.identifier
            case 'VISA_Instrument:ErrorChecking'
                % Perform your own additional steps here
                rethrow(ME);
            otherwise
                rethrow(ME);
        end
    end
    if(sum(OFDMSymbolNumArray)~=0)
        print(num2str(sum(OFDMSymbolNumArray)));
        figure(3)
        scatter(real(OFDMSymbolWithoutPaddingAcc(1,:)),imag(OFDMSymbolWithoutPaddingAcc(1,:)),'.');
        xlim([-1.5,1.5]);
        ylim([-1.5,1.5]);
        figure(4)
        scatter(real(OFDMSymbolWithoutPaddingAcc(20,:)),imag(OFDMSymbolWithoutPaddingAcc(20,:)),'.');
        xlim([-1.5,1.5]);
        ylim([-1.5,1.5]);
        figure(5)
        scatter(real(OFDMSymbolWithoutPaddingAcc(40,:)),imag(OFDMSymbolWithoutPaddingAcc(40,:)),'.');
        xlim([-1.5,1.5]);
        ylim([-1.5,1.5]);
        figure(6)
        scatter(real(OFDMSymbolWithoutPaddingAcc(58,:)),imag(OFDMSymbolWithoutPaddingAcc(58,:)),'.');
        xlim([-1.5,1.5]);
        ylim([-1.5,1.5]);
    end
	if (sum(OFDMSymbolNumArray) >= TotalFrameNum)
	    disp('Receiving Done, Begin to calculate ....')
	    break;
	end
end
% Drawing the receiveing constellation
scatterplot(OFDMSymbolWithoutPaddingAcc(SpecificSubcarrier,:))
%% SNR calculation
% First part: Ideal normalization factor calculation
% normalization scaling factor for ideal symbol calculation 
DataModIdeal = [];
qamMapper_Ideal.ModulationOrder = 2^M;
for nSymbol = 1:2^M
	DataModIdeal = [DataModIdeal, qammod(nSymbol-1,2^M,'gray','InputType','integer','UnitAveragePower',true)];
end
IdealScaleFactor = sqrt((2^M)/(sum(abs(DataModIdeal).^2)));
NormalizedDataModIdeal = DataModIdeal./IdealScaleFactor;
AverageEnergyIdeal = sum((abs(NormalizedDataModIdeal)).^2)/(2^M);
%%% % The bit sequence generation
%%% [DataLogic,~] = prbs_gen(POLY_LENGTH, POLY_TAP, NBITS, NUM);
%%% DataSerial = double(DataLogic);
%%% % Symbol Mapping
%%% DataSymbol =  qammod(DataSerial',2^M,'gray','InputType','bit','UnitAveragePower',true);
%%% QuantizerInst = quantizer('fixed','Nearest','saturate',[WORD_LENGTH,FRACTION_LENGTH]);
%%% DataSymbolQuan = quantize(QuantizerInst,DataSymbol);
%%% DataSymbolQuanReshape = reshape(DataSymbolQuan,SubCarrierNum,FrameNum);
ModulationDataAll = [];
for nReceive = 1:length(OFDMSymbolNumArray)
	ModulationDataAll = [ModulationDataAll,DataSymbolQuanReshape(:,1:OFDMSymbolNumArray(nReceive))];
end

% Source normalization coefficient calculation
for nSubCarrier = 1:SubCarrierNum
	SourceScaleFactor(nSubCarrier) = sqrt(sum(OFDMSymbolNumArray)/sum(abs(ModulationDataAll(nSubCarrier,:)).^2));
end
% Received normalization coefficient calculation
for nSubCarrier = 1:SubCarrierNum
	ReceivedScaleFactor(nSubCarrier) = sqrt(sum(OFDMSymbolNumArray)/sum(abs(OFDMSymbolWithoutPaddingAcc(nSubCarrier,:)).^2));
end

% Normalization source and destination information
for nSubCarrier = 1:SubCarrierNum
	NormalizedSourceSymbol = ModulationDataAll(nSubCarrier,:).*SourceScaleFactor(nSubCarrier);
	NormalizedReceivedSymbol = OFDMSymbolWithoutPaddingAcc(nSubCarrier,:).*ReceivedScaleFactor(nSubCarrier);
	NoiseEnergy(nSubCarrier) = sum(abs(NormalizedReceivedSymbol-NormalizedSourceSymbol).^2)/sum(OFDMSymbolNumArray);
	EvmResult(nSubCarrier) = sqrt(NoiseEnergy(nSubCarrier)/AverageEnergyIdeal);
end
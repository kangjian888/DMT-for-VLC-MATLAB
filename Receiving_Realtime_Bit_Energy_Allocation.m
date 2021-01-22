% Starting data acquisition and BER calculation
%% Global Parameter
N= 128;
Ncp = 32;
RepeatTime = 2;
FrameNum = 200; % ofdm symbol number of each package transmitted
TotalFrameNum = 2000; % total ofdm symbol number you want to get
SubCarrierNum = 59; % every OFDM frame has 52 subcarriers
LowPaddingNum = 0; % padding number of low frequency
HighPaddingNum = 0; % padding number of high frequency
Pilot_pos = [7,21,43,57]; % this is the position of pilot
ThresholdValue = 60;

POLY_LENGTH = 7;
POLY_TAP = 1;
PILOT_NUM = 4; % there is four subcarriers are used for pilot

WORD_LENGTH = 14;
FRACTION_LENGTH = 12;
% RArray = 4.*ones(1,59);
% SArray = 1.*ones(1,59);
%RArray = 6.*[zeros(1,3),ones(1,47),zeros(1,9)];
%SArray = 1.*[zeros(1,3),ones(1,47),zeros(1,9)];
RSum = 236;
% Data for 30cm, 10 LED result
RSum = 94
;
if RSum ==59
    RArray = [0,0,0,0,0,1,2,2,2,2,2,2,3,3,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
    SArray = [0,0,0,0,1.34740763690093,1.53810971761448,2.0,1.61166918092065,1.39550531051351,1.37639393769851,1.15752738516162,1.28422146294437,2.01177654843575,1.95742022281563,1.10943813274643,2.0,1.31772081388189,1.26525780443209,1.15553032362836,1.37398185185999,1.32114978600997,1.27756346748179,1.53480679215627,1.94283247611855,1.61198287167926,1.70486791912331,1.80944055492069,1.90864429728138,2.0,1.09654258955452,1.21967494460627,1.24846503544627,1.29813816636718,1.25336225441653,1.49933569758692,1.52325678492961,1.69714189584405,1.87525874883677,1.02808084782147,1.06773242907425,1.14568365321556,1.29895854431691,1.39216647784539,1.27614270513430,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
end
if RSum ==94
     RArray = [0,0,0,0,1,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,3,3,3,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0];
     SArray = [0,0,0,1.20695915231454,1.30829373132451,1.49345992002292,0.973179960265440,1.56488402516189,1.35499517723836,1.33643858859650,1.12392551472139,1.24694178924714,0.976688336574151,0.950299128839590,1.07723224554594,0.984960621637173,1.27946868729546,1.22852862690850,1.12198642590113,1.33409652321441,1.28279811976801,1.24047706878073,1.49025287521060,0.943216989390745,1.56518860978931,1.65537729655104,1.75691429258327,0.926619124411682,0.987208038183439,1.06471104716620,1.18426899223466,1.21222333532928,1.26045450459040,1.21697839289631,1.45581147128993,1.47903815322332,1.64787555205480,0.910410924847305,0.998236681805679,1.03673721703304,1.11242559455160,1.26125106778297,1.35175326756677,1.23909754975322,1.62901610245745,1.69335970789085,0.904654300341347,1.13929390797981,1.79910582981138,0.992905529915808,0,0,0,0,0,0,0,0,0];
end
if RSum == 118
    RArray = [0,0,0,1,1,2,3,3,3,3,4,3,4,4,4,4,3,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,2,3,2,2,2,2,2,2,2,1,1,2,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0]; %59 used subcarrier
    SArray = [0,0,0,1.60114538318393,0.867787639598932,0.990607864122556,1.29101519086046,1.03798327694079,0.898764580438056,0.886456046145664,1.49099341554999,0.827093064882725,1.29566938360459,1.26066160554595,1.42905038108800,1.30664335154016,0.848668067043935,1.62975932972685,1.48842103096996,0.884902560608703,0.850876470460829,0.822805033548557,0.988480640104630,1.25126647824621,1.03818530708781,1.09800721537156,1.16535642603357,1.22924784171574,1.30962476198826,1.41243982808031,1.57104474142988,1.60812881946189,0.836056011865913,1.61443685271689,0.965635751455249,0.981041945804497,1.09303132884592,1.20774613318446,1.32425530005135,1.37532989864998,1.47573768460704,0.836584370044297,0.896614230654711,1.64378000473114,1.08052190771337,1.12320084451895,1.20010942672860,1.51138104164776,1.19334195682501,1.31718302323142,0,0.872102222461623,1.06862719141565,1.00219713746568,0,0,0,0,0];
end
if RSum == 177
    RArray = [0,0,0,2,3,3,4,4,4,4,5,5,5,5,5,5,5,5,5,4,5,5,4,4,4,4,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,2,3,2,2,2,2,2,2,1,2,1,1,0,0,0,0,0];
    SArray = [0,0,0.992458235232082,1.32788791657893,1.43937550316453,0.821547016686675,1.07068570416414,0.860837164155468,0.745378051704207,0.735170137963630,1.23653489620695,1.37187653073704,1.07454559494686,1.04551237537882,1.18516329195278,1.08364670442237,1.40766239374601,1.35161849985952,1.23440152441168,0.733881776085454,1.41132541178676,1.36476409103687,0.819782832685069,1.03772065568423,0.861004715081095,0.910617192493179,0.966472334673574,1.01945980211314,1.08611929619947,1.17138755811982,1.30292436305889,1.33367953347212,1.38674312451536,1.33891101415130,0.800836738275590,0.813613653808488,0.906490509394283,1.00162765567325,1.09825293185246,1.14061094820471,1.22388276543202,1.38761949770179,0.743594688666901,1.36324635396697,0.896115992913912,0.931511182552388,0.995294258136900,1.25344309369313,0.989681749987943,1.09238763627686,1.05331688725751,1.44653198315519,0.886251273439487,0.831158421244866,0.822659667861146,0.800948722764478,0.791828145271692,0,0];
end
if RSum == 236
    RArray = [0,0,1,3,4,4,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6,5,5,5,5,5,5,5,5,5,5,5,5,4,4,4,4,4,4,4,4,4,4,3,3,3,3,3,3,2,3,2,2,1,1,1,1,0];
    SArray = [0,0.908726341594383,0.913292738626260,1.22196617335489,1.32456072057207,0.756014609123124,0.985280230694427,0.792170696209663,0.685921420152911,1.35305552379247,1.13790011677388,1.26244594415697,0.988832229255790,0.962114904869599,1.09062627544599,0.997207369793046,1.29537727332663,1.24380383727276,1.13593691781976,1.35068433776914,1.29874810311000,1.25590084301167,0.754391148929969,0.954944707912333,0.792324882087005,0.837979916986403,0.889379657506833,0.938140469330246,0.999482730141165,1.07794939169184,1.19899388963472,1.22729580987326,1.27612667314268,1.23210999060848,0.736956329171072,0.748714067445578,0.834182407353212,0.921730741162832,1.01064850109416,1.04962774211791,1.12625720954142,1.27693314059802,1.36856062155951,1.25450417139781,0.824635435771012,0.857207254444632,0.915902540261835,1.15345959669877,0.910737725505953,1.00525106302681,0.969296873615786,1.33114634904243,0.815557596175605,0.764859339993781,0.757038507361048,0.737059380984145,0.728666325336528,1.01852458882339,0.742856616946392];
end
if RSum == 295
    RArray = [0,1,2,4,5,5,6,6,7,7,7,7,7,7,7,7,7,7,7,7,7,7,6,6,6,6,6,6,6,6,6,6,6,6,5,5,5,5,5,5,5,5,5,5,4,4,4,4,4,4,3,4,3,3,2,2,2,2,1];
    SArray = [0.775680328768165,0.886473254458217,0.890927828572240,1.19204240149041,1.29212459125252,0.737501159910775,0.961152475369633,0.772771848919133,1.33824885621498,1.31992160757057,1.11003497267957,1.23153089491263,0.964617491822218,0.938554427052369,1.06391878339959,0.972787540107140,1.26365579457323,1.21334529997256,1.10811984896592,1.31760848766327,1.26694407874802,1.22514607161934,0.735917455361852,0.931559815428520,0.772922259059374,0.817459283592704,0.867600336153709,0.915167082674214,0.975007180939080,1.05155233391578,1.16963266802101,1.19724152471735,1.24487660724891,1.20193781475382,0.718909583238936,0.730379395476340,0.813754768196502,0.899159199478246,0.985899522078732,1.02392423102028,1.09867717947683,1.24566332575483,1.33504701315100,1.22378360200195,0.804441584918636,0.836215777853955,0.893473720821076,1.12521342869582,0.888435383165663,0.980634257641800,0.945560522195027,1.29854894942897,0.795586045422036,0.746129299099480,0.738499984707264,0.719010111496930,0.710822586810743,0.993582738468699,0.724665383491465];
end
if RSum == 354
    RArray = [1,2,3,5,6,6,7,7,8,8,8,8,8,8,8,8,8,8,8,8,8,8,7,7,7,7,7,7,7,7,7,7,7,7,6,6,6,6,6,6,6,6,6,6,5,5,5,5,5,5,4,5,4,4,3,3,3,3,2];
    SArray = [0.775680328768165,0.886473254458217,0.890927828572240,1.19204240149041,1.29212459125252,0.737501159910775,0.961152475369633,0.772771848919133,1.33824885621498,1.31992160757057,1.11003497267957,1.23153089491263,0.964617491822218,0.938554427052369,1.06391878339959,0.972787540107140,1.26365579457323,1.21334529997256,1.10811984896592,1.31760848766327,1.26694407874802,1.22514607161934,0.735917455361852,0.931559815428520,0.772922259059374,0.817459283592704,0.867600336153709,0.915167082674214,0.975007180939080,1.05155233391578,1.16963266802101,1.19724152471735,1.24487660724891,1.20193781475382,0.718909583238936,0.730379395476340,0.813754768196502,0.899159199478246,0.985899522078732,1.02392423102028,1.09867717947683,1.24566332575483,1.33504701315100,1.22378360200195,0.804441584918636,0.836215777853955,0.893473720821076,1.12521342869582,0.888435383165663,0.980634257641800,0.945560522195027,1.29854894942897,0.795586045422036,0.746129299099480,0.738499984707264,0.719010111496930,0.710822586810743,0.993582738468699,0.724665383491465];
end
%% Synchronization
i = 1; % index of input data
% Some registers used in Frame Syn
L1 = 32;
L2 = 70;
L3 = N;
L4 = N*(RepeatTime-1) + Ncp;
SynchronizationOffset = 0;% this parameter is the offset position of the beginning of the sequence. 
                                  % if this value is 0, Then the sequence beginning posision is the first sample after the repeaded time peak
SpecificSubcarrier = 40;
% Parameter Check
if SynchronizationOffset > Ncp
    error('The SynchronizationOffset is invalid, Please adjust the parameters.');
end

% globle registers and varibles
OFDMSymbolWithoutPaddingAcc = [];
OFDMSymbolNumArray = [];
ErrorBitNumArray = [];
ErrorBitRatioArry = [];
% Getting Data Processing
[~,~,DataSerialOri]= TransmitterTopNew_Bit_Energy_Allocation(N,RArray,RSum,SArray,Ncp,RepeatTime,FrameNum,LowPaddingNum,HighPaddingNum,SubCarrierNum,Pilot_pos,WORD_LENGTH,FRACTION_LENGTH);
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
        waveformASC_NEG = scope.QueryASCII_ListOfDoubles('FORM ASC;:CHAN1:DATA?', samplesCount); % samplesCount is the maximum allowed samples to read
        waveformASC_POS = scope.QueryASCII_ListOfDoubles('FORM ASC;:CHAN2:DATA?', samplesCount); % samplesCount is the maximum allowed samples to read
        waveformASC_TRAN = scope.QueryASCII_ListOfDoubles('FORM ASC;:CHAN3:DATA?', samplesCount); % samplesCount is the maximum allowed samples to read
        toc
        fprintf('Samples count: %d\n', size(waveformASC_POS, 2));
        fprintf('Samples count: %d\n', size(waveformASC_NEG, 2));
        scope.ErrorChecking(); % Error Checking after the data transfer
        scope.Write('CHAN2:STAT ON'); % Switch Channel 1 ON
        scope.Write('CHAN1:STAT ON'); % Switch Channel 1 ON
        scope.Write('CHAN3:STAT ON'); % Switch Channel 1 ON
        waveformASC = waveformASC_POS-waveformASC_NEG;
        %waveformASC = waveformASC_TRAN;
        %waveformASC = -waveformASC_POS;
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
                if BitSumPro > 50 % first threshold
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
                    PayLoadLength = LengthReceivedSeq + L1 + L2 - i + 2;
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
                            OFDMSymbolNumArray = [OFDMSymbolNumArray,FrameCounter-RepeatTime]; % the OFDM number getted each sampling time
                            % demodulation and calculate the bit error rate
                            DataSerialReshape = [];
                            for nSubcarrier = 1:SubCarrierNum
                                if RArray(nSubcarrier) == 0
                                    DataSerialReshape = DataSerialReshape;
                                else
                                    OFDMSymbolWithoutPaddingOneSubcarrier = OFDMSymbolWithoutPadding(nSubcarrier,:);
                                    DataSerialOneSubcarrier =  qamdemod(OFDMSymbolWithoutPaddingOneSubcarrier ./ SArray(nSubcarrier),2^RArray(nSubcarrier),'gray','OutputType','bit','UnitAveragePower',true);
                                    DataSerialReshape = [DataSerialReshape;DataSerialOneSubcarrier];
                                end
                            end
                            DataSerial = reshape(DataSerialReshape,1,RSum*(FrameCounter-RepeatTime));
                            [ErrorNum, ErrorRate] = biterr(DataSerialOri(:,1:RSum*(FrameCounter-RepeatTime)),DataSerial);
                            ErrorBitNumArray = [ErrorBitNumArray,ErrorNum];
                            ErrorBitRatioArry = [ErrorBitRatioArry,ErrorRate];
                            disp('Bit Error Rate in this time is');
                            disp(ErrorRate);
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
	% Drawing the receiveing constellation
    if(sum(OFDMSymbolNumArray)~=0)
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
        scatter(real(OFDMSymbolWithoutPaddingAcc(59,:)),imag(OFDMSymbolWithoutPaddingAcc(59,:)),'.');
        xlim([-1.5,1.5]);
        ylim([-1.5,1.5]);
    end
	if (sum(OFDMSymbolNumArray) >= TotalFrameNum)
	    disp('Receiving Done, Begin to calculate ....')
        % calculate the final bit error rate
        TotalBitReceived = sum(OFDMSymbolNumArray)*RSum;
        TotalErrorBitReceived = sum(ErrorBitNumArray);
        BitErrorRate = TotalErrorBitReceived/TotalBitReceived;
        disp('Bit Error Rate is')
        disp(BitErrorRate);
	    break;
	end
end


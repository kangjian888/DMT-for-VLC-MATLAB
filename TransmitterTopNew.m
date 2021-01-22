function [Output2DAC,DataSymbolQuanReshape,DataSerial]= TransmitterTopNew(N,M,Ncp,RepeatTime,FrameNum,LowPaddingNum,HighPaddingNum,SubCarrierNum,Pilot_pos,WORD_LENGTH,FRACTION_LENGTH)
%clear all;
%close all;
%clc;

%%% % Parameters used for test
%%% N = 128;
%%% M = 4;
%%% Ncp = 32;
%%% RepeatTime = 2;
%%% FrameNum = 100;
%%% LowPaddingNum = 0;
%%% HighPaddingNum = 0;
%%% SumcarrierNum = 59;
%%% Pilot_pos = [7,21,43,57];
%%% WORD_LENGTH = 14;
%%% FRACTION_LENGTH = 12;

%% Parameter Aera
POLY_LENGTH = 7;
POLY_TAP = 1;
PILOT_NUM = 4; % there is four subcarriers are used for pilot
%Pilot_pos = [7,21,43,57]; % this is the position of pilot
NBITS = 1; %Data width of parallel output
NUM = FrameNum*SubCarrierNum*M; %Byte number

%% Parameters Checking
if N ~= (SubCarrierNum + LowPaddingNum + HighPaddingNum + PILOT_NUM + 1) * 2
	error('The number of subcarrier is invalid');
else
	disp('Parameter checed, simulation begins ...');
end

[DataLogic,~] = prbs_gen(POLY_LENGTH, POLY_TAP, NBITS, NUM);
DataSerial = double(DataLogic);

%% Symbol Mapping
DataSymbol =  qammod(DataSerial',2^M,'gray','InputType','bit','UnitAveragePower',true);
DataSymbolReal = real(DataSymbol);
DataSymbolImag = imag(DataSymbol);
[DataSymbolRealQuan,~,DataSymbolRealBinary] = DACInputGen(DataSymbolReal, WORD_LENGTH, FRACTION_LENGTH);
[DataSymbolImagQuan,~,DataSymbolImagBinary] = DACInputGen(DataSymbolImag, WORD_LENGTH, FRACTION_LENGTH);
QuantizerInst = quantizer('fixed','Nearest','saturate',[WORD_LENGTH,FRACTION_LENGTH]);
DataSymbolQuan = quantize(QuantizerInst,DataSymbol);

DataSymbolQuanReshape = reshape(DataSymbolQuan,SubCarrierNum,FrameNum);
DataSymbolAddingCP = [];

% register used in scrambler
ScramblerRegister = ones(1,7);

for nFrame = 1:FrameNum
	% scrambler data generation
	ScramblerData = xor(ScramblerRegister(4),ScramblerRegister(7));
	ScramblerRegister = [ScramblerData,ScramblerRegister(2:7)];
	% pilot generation
	if ScramblerData == 0
		PilotSeq = [1,-1,1,1];
    elseif ScramblerData ==1
		PilotSeq = [-1,1,-1,-1];
	else
		error('ScramblerData is invalid');
    end
	% Pick up one symbol to process
	DataSymbolOneFrame =  DataSymbolQuanReshape(:,nFrame);
	% Hermitian Symmetry and IFFT transformation
	% Adding Padding and pilot to the symbol
	DataSymbolOneFrameAddPadiing = [zeros(LowPaddingNum,1);DataSymbolOneFrame;zeros(HighPaddingNum,1)];
	DataSymbolOneFrameAddPilot = [DataSymbolOneFrameAddPadiing(1-0:Pilot_pos(1)-1-0);PilotSeq(1);DataSymbolOneFrameAddPadiing(Pilot_pos(1)+1-1:Pilot_pos(2)-1-1);PilotSeq(2);DataSymbolOneFrameAddPadiing(Pilot_pos(2)+1-2:Pilot_pos(3)-1-2);PilotSeq(3);DataSymbolOneFrameAddPadiing(Pilot_pos(3)+1-3:Pilot_pos(4)-1-3);PilotSeq(4);DataSymbolOneFrameAddPadiing(Pilot_pos(4)+1-4:end)];
	DataSymbolOneFrameHermitianSymmetry = [0;DataSymbolOneFrameAddPilot;0;flipud(conj(DataSymbolOneFrameAddPilot))];

	%% IFFT transformation
	DataSymbolOneFrameAfterIFFT = ifft(DataSymbolOneFrameHermitianSymmetry,N,1); %Every coloum one IFFT symbol

	%% Adding CP
	DataSymbolOneFrameAddingCP = [DataSymbolOneFrameAfterIFFT(N-Ncp+1:N);DataSymbolOneFrameAfterIFFT];
	DataSymbolAddingCP = [DataSymbolAddingCP,DataSymbolOneFrameAddingCP];
    end
QuantizerInst = quantizer('fixed','Nearest','saturate',[WORD_LENGTH,FRACTION_LENGTH]);
DataSymbolAddingCPQuan = quantize(QuantizerInst,reshape(DataSymbolAddingCP,1,size(DataSymbolAddingCP,1)*size(DataSymbolAddingCP,2)));
%% Short Training Sequence Generation
[ShortTrainingSeqOutput] = ShortTrainingGen();
%% Long Training Sequence Generation
[LongTrainingSeqOutput,~,~] = TrainingSeqGenParkMethod(N,Ncp,RepeatTime);
%% Adding Training Sequence to output signal
Output2DAC = [ShortTrainingSeqOutput,LongTrainingSeqOutput,reshape(DataSymbolAddingCP,1,size(DataSymbolAddingCP,1)*size(DataSymbolAddingCP,2))];
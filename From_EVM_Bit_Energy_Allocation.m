clc;
clear all;
close all;
path = 'D:\OneDrive - HKUST Connect\Lab\OFDM_FPGA_VLC_Implementation\measurement\59subcarrier_measure\Rece\20Mbps_2bit_SNR_10LED_30cm\';
% read EVM result
T = readtable([path,'EVM.xlsx'],'Sheet','Sheet1');
EvmResult = table2array(T(2,2:60));

% Parameters
%TargetSpeed = 4e7;
AverageEnergyIdeal = 1;
SAMPLING_RATE = 2e7;
N= 128;
Ncp = 32;
SubCarrierNum = 59;
% SNR linear and dB calculation
%R_T = ceil(TargetSpeed*(1/SAMPLING_RATE*(N+Ncp)));
R_T = 94;
SNR = 1./EvmResult.^2;
SNR_dB = 10*log10(SNR);
NoiseVarance = AverageEnergyIdeal./SNR;
D = SubCarrierNum; % number of used channel, right side is initial value
L = 1:D; % index of used channel, right side is initial value
R_MAX = 6; % Maxmum bit in one OFDM Frame
R = zeros(1,SubCarrierNum); % bit number on each subcarrier before quantification
RQ = zeros(1,SubCarrierNum); % bit number on each subcarrier after quantification
RDelta = zeros(1,SubCarrierNum); % quantification error
LDN = log2(NoiseVarance);
Flag = 1;
while(Flag == 1)
	Flag = 0;
	for i = 1:D
		R(L(i)) = (R_T + sum(LDN(L)))/D - LDN(L(i));
		if(R(L(i)) <= 0)
			R(L(i)) =0;
			D = D-1;
			L(i) = [];
			Flag = 1;
			break; % break out this turn to next turn
		end
	end
end

for i = 1:D
	if (R(L(i)) >= R_MAX-0.5)
		RQ(L(i)) = R_MAX;
	elseif (R(L(i))>= 0.5 && R(L(i)) < R_MAX-0.5)
		RQ(L(i)) = round(R(L(i)));
	else
		RQ(L(i)) = 0;
	end
	RDelta(L(i)) = R(L(i)) - RQ(L(i));
end

R_SUM = sum(RQ);
while(R_SUM ~= R_T)
	if (R_SUM > R_T)
		[MIN_VALUE,MIN_INDEX] = min(RDelta);
		RQ(MIN_INDEX) = RQ(MIN_INDEX) - 1;
		R_SUM = R_SUM -1;
		RDelta(MIN_INDEX) = RDelta(MIN_INDEX) + 1;
	elseif (R_SUM < R_T)
		[MAX_VALUE,MAX_INDEX] = max(RDelta);
		RQ(MAX_INDEX) = RQ(MAX_INDEX) + 1;
		R_SUM = R_SUM + 1;
		RDelta(MAX_INDEX) = RDelta(MAX_INDEX) - 1;	
	end
end


% Energy distribution
S_T = AverageEnergyIdeal * (SubCarrierNum);
S = zeros(1,SubCarrierNum);
for i = 1:D
	S(L(i)) = (S_T*NoiseVarance(L(i))*(2^RQ(L(i))))/(sum(NoiseVarance(L).*2.^RQ(L)));
end

% Output 
S_Array = S;
R_Array = RQ;

%Speed Calculation
DataRate = (sum(R_Array))/(1/SAMPLING_RATE*(N+Ncp));
%% plot and save the result
% generate coe file used by FPGA
fid_w_1 = fopen([path,num2str(DataRate/1e6),'Mpbs','Bit_Allocation','.coe'],'w');
fprintf(fid_w_1,'%s\n','memory_initialization_radix=10;');
fprintf(fid_w_1,'%s','memory_initialization_vector=');
[row,col] = size(R_Array);
for i = 1:col
    if i == col
        fprintf(fid_w_1,'%g',R_Array(1,i));
        fprintf(fid_w_1,'%s',';');
    else
        fprintf(fid_w_1,'%g',R_Array(1,i));
        fprintf(fid_w_1,'%s',',');
    end
end
fclose(fid_w_1);
% transfer energy into fixed-number for FPGA
WORD_LENGTH = 14;
FRACTION_LENGTH = 12;
[~,~,EnergyBinary] = DACInputGen(S_Array, WORD_LENGTH, FRACTION_LENGTH);
fid_w_2 = fopen([path,num2str(DataRate/1e6),'Mpbs','Energy_Allocation','.coe'],'w');
fprintf(fid_w_2,'%s\n','memory_initialization_radix=2;');
fprintf(fid_w_2,'%s','memory_initialization_vector=');
[row,col] = size(EnergyBinary);
for i = 1:row
	for j = 1:col
		if(j == col)
			fprintf(fid_w_2,'%g',EnergyBinary(i,j));
            if(i == row)
                fprintf(fid_w_2,'%s',';');
            else
                fprintf(fid_w_2,'%s',',');
            end
		else
			fprintf(fid_w_2,'%g',EnergyBinary(i,j));
		end
	end
end
fclose(fid_w_2);
name = ["Subcarrier Index";"EVM";"SNR (linear)";"SNR (dB)";"Bit Number";"Energy Allcation"];
data = [1:SubCarrierNum;EvmResult;SNR;SNR_dB;R_Array;S_Array];
result_2 = table(name,data);
writetable(result_2,[path,'Result with',num2str(DataRate/1e6),'Mpbs Daterate.xlsx']);  

%% Result Plot
figure(1)
plot(1:SubCarrierNum,SNR,'LineWidth',2);
%legend({'Linear SNR'},'Location','southeast','FontSize',10);
title('SNR for different subcarrier');
xlabel('Index of Subcarrier');
ylabel('SNR');
set(gca, 'fontsize', 16);
set(gca, 'XMinorTick', 'on');
set(gca, 'YMinorTick', 'on');
set(gca, 'XGrid', 'on');
set(gca, 'YGrid', 'on');
set(gca, 'LineWidth', 1.5);
saveas(1,fullfile(path,'SNR(linear)'),'tif');

figure(2)
plot(1:SubCarrierNum,SNR_dB,'LineWidth',2);
title('SNR(dB) for different subcarrie');
xlabel('Index of Subcarrier');
ylabel('SNR(dB)');
set(gca, 'fontsize', 16);
set(gca, 'XMinorTick', 'on');
set(gca, 'YMinorTick', 'on');
set(gca, 'XGrid', 'on');
set(gca, 'YGrid', 'on');
set(gca, 'LineWidth', 1.5);
saveas(2,fullfile(path,'SNR(dB)'),'tif');

figure(3)
plot(1:SubCarrierNum,RQ,'LineWidth',2);
%legend({'Bit Number on Each Subcarrier'},'Location','southeast','FontSize',10);
title('Number of Bit on Each Subcarrier');
xlabel('Index of Subcarrier');
ylabel('Number of Bit');
set(gca, 'fontsize', 16);
set(gca, 'XMinorTick', 'on');
set(gca, 'YMinorTick', 'on');
set(gca, 'XGrid', 'on');
set(gca, 'YGrid', 'on');
set(gca, 'LineWidth', 1.5);
saveas(3,fullfile(path,['Bit Number_2',num2str(DataRate/1e6),'Mpbs']),'tif');

figure(4)
plot(1:SubCarrierNum,S,'LineWidth',2);
%legend({'Allocated Energy on Each Subcarrier'},'Location','southeast','FontSize',10);
title('Energy Allocation Scaling Factor on Each Subcarrier');
xlabel('Index of Subcarrier');
ylabel('Energy Scaling Factor');
set(gca, 'fontsize', 16);
set(gca, 'XMinorTick', 'on');
set(gca, 'YMinorTick', 'on');
set(gca, 'XGrid', 'on');
set(gca, 'YGrid', 'on');
set(gca, 'LineWidth', 1.5);
saveas(4,fullfile(path,['Energy Scaling Factor_2',num2str(DataRate/1e6),'Mpbs']),'tif');

figure(5)
scatter(1:SubCarrierNum,RQ,'LineWidth',2);
xlim([0,60]);
ylim([0,8]);
for subcarrierindex = 1:SubCarrierNum
    line([subcarrierindex,subcarrierindex],[RQ(subcarrierindex),0]);
end
%legend({'Bit Number on Each Subcarrier'},'Location','southeast','FontSize',10);
title('Number of Bit on Each Subcarrier');
xlabel('Index of Subcarrier');
ylabel('Number of Bit');
set(gca, 'fontsize', 16);
set(gca, 'XMinorTick', 'on');
set(gca, 'YMinorTick', 'on');
set(gca, 'XGrid', 'on');
set(gca, 'YGrid', 'on');
set(gca, 'LineWidth', 1.5);
saveas(5,fullfile(path,['Bit Number_3',num2str(DataRate/1e6),'Mpbs']),'tif');

figure(6)
scatter(1:SubCarrierNum,S,'LineWidth',2);
xlim([0,60]);
ylim([0,2]);
for subcarrierindex = 1:SubCarrierNum
    line([subcarrierindex,subcarrierindex],[S(subcarrierindex),0]);
end
%legend({'Allocated Energy on Each Subcarrier'},'Location','southeast','FontSize',10);
title('Energy Allocation Scaling Factor on Each Subcarrier');
xlabel('Index of Subcarrier');
ylabel('Energy Scaling Factor');
set(gca, 'fontsize', 16);
set(gca, 'XMinorTick', 'on');
set(gca, 'YMinorTick', 'on');
set(gca, 'XGrid', 'on');
set(gca, 'YGrid', 'on');
set(gca, 'LineWidth', 1.5);
saveas(6,fullfile(path,['Energy Scaling Factor_3',num2str(DataRate/1e6),'Mpbs']),'tif');
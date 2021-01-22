clc;
clear all;
close all;
path = 'D:\OneDrive - HKUST Connect\Lab\OFDM_FPGA_VLC_Implementation\measurement\59subcarrier_measure\Rece\20Mbps_2bit_SNR_10LED_30cm\';
SubCarrierNum = 59;
BitOneEachSubcarrier = 2;
EnergyScalingFactor = 1;
figure(1)
scatter(1:SubCarrierNum,BitOneEachSubcarrier*ones(1,SubCarrierNum),'LineWidth',2);
xlim([0,60]);
ylim([0,8]);
for subcarrierindex = 1:SubCarrierNum
    line([subcarrierindex,subcarrierindex],[BitOneEachSubcarrier,0]);
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
saveas(1,fullfile(path,['Bit Number_3',num2str(BitOneEachSubcarrier),'bit']),'tif');

figure(2)
scatter(1:SubCarrierNum,EnergyScalingFactor*ones(1,SubCarrierNum),'LineWidth',2);
xlim([0,60]);
ylim([0,2]);
for subcarrierindex = 1:SubCarrierNum
    line([subcarrierindex,subcarrierindex],[EnergyScalingFactor,0]);
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
saveas(2,fullfile(path,['Energy Scaling Factor',num2str(EnergyScalingFactor),'factor']),'tif');
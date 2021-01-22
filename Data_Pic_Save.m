close all;
path =  'D:\OneDrive - HKUST Connect\Lab\OFDM_FPGA_VLC_Implementation\measurement\59subcarrier_measure\Rece\20Mbps_3bit_SNR_10LED_30cm\';
for i = 1:59
        figure(i)
        scatter(real(OFDMSymbolWithoutPaddingAcc(i,:)),imag(OFDMSymbolWithoutPaddingAcc(i,:)),'.');
        xlim([-1.5,1.5]);
        ylim([-1.5,1.5]);
        box on;
        saveas(i,[path,num2str(i),'.jpg']);
 end
title = ["Subcarrier Index";"EVM"];
data = [1:SubCarrierNum;EvmResult];
result =table(title,data);
writetable(result,[path,'EVM.xlsx']);   
close all;
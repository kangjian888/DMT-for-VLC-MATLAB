clear;
close all;
clc;
try
    scope = VISA_Instrument('TCPIP::192.168.168.1::INSTR'); % Adjust the VISA Resource string to fit your instrument
    scope.SetTimeoutMilliseconds(3000); % Timeout for VISA Read Operations
    % scope.AddLFtoWriteEnd = false;
catch ME
    error ('Error initializing the instrument:\n%s', ME.message);
end

% Setting the oscilloscope, including time span, sampling
% rate, channel setting, timeout setting.
try
    idnResponse = scope.QueryString('*IDN?');
    fprintf('\nInstrument Identification string: %s\n', idnResponse);
    %scope.Write('*RST;*CLS'); % Reset the instrument, clear the Error queue
    scope.Write('SYST:DISP:UPD ON'); % Display update ON - switch OFF after debugging
    scope.ErrorChecking(); % Error Checking after Initialization block
    %-----------------------------------------------------------
    % Basic Settings:
    %-----------------------------------------------------------
    scope.Write('ACQ:POIN:AUTO RECL'); % Define Horizontal scale by number of points
    scope.Write('TIM:RANG %f',2e-3); % 100ms Acquisition time
    scope.Write('ACQ:POIN %d', 40000); % 100ksamples
    %scope.Write('CHAN3:POS 0'); % Offset 0
    %scope.Write('CHAN3:STAT ON'); % Switch Channel 1 ON
    scope.Write('CHAN2:POS 0'); % Offset 0
    scope.Write('CHAN2:STAT ON'); % Switch Channel 1 ON
    scope.Write('CHAN1:POS 0'); % Offset 0
    scope.Write('CHAN1:STAT ON'); % Switch Channel 1 ON
    scope.ErrorChecking(); % Error Checking after Basic Settings block

    %-----------------------------------------------------------
    % Trigger Settings:
    %-----------------------------------------------------------
    %scope.Write('TRIG1:MODE AUTO'); % Trigger Auto mode in case of no signal is applied
    %scope.Write('TRIG1:SOUR CHAN1'); % Trigger source CH1
    %scope.Write('TRIG1:TYPE EDGE;:TRIG1:EDGE:SLOP POS'); % Trigger type Edge Positive
    %scope.Write('TRIG1:LEV1 0.04'); % Trigger level 40mV
    %scope.QueryString('*OPC?'); % Using *OPC? query waits until all the instrument settings are finished
    %scope.ErrorChecking(); % Error Checking after Trigger Settings block
    
catch ME
    switch ME.identifier
        case 'VISA_Instrument:ErrorChecking'
            % Perform your own additional steps here
            rethrow(ME);
        otherwise
            rethrow(ME);
    end
end 
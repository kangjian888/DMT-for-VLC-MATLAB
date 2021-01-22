clc;
clear all;
close;
try
    scope.Close(); % Adjust the VISA Resource string to fit your instrument
    scope.SetTimeoutMilliseconds(3000); % Timeout for VISA Read Operations
    % scope.AddLFtoWriteEnd = false;
catch ME
    error ('Error initializing the instrument:\n%s', ME.message);
end
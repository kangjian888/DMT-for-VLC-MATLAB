classdef VISA_Instrument < handle
% VISA_Instrument This class defines instrument remote-control connection over VISA.NET. Version 1.7.0
% MATLAB Instrument Control Toolbox is not required

% Version 1.8.0
% - Added varargin for passing cmd/query sprintf parameters to the following functions:
% 	WriteWithSRQsync
% 	QueryWithSRQsync
% 	WriteWithSTBpollSync
% 	QueryWithSTBpollSync
% 	QueryBinaryDataBlock
% 	QueryASCII_ListOfDoubles

% Version 1.7.0
% - Session.Clear() in constructor is not performed on SERIAL a SOCKET session.
% - Added settings for SERIAL session type - termination character enable in constructor
% - ReadBinaryDataBlock() - reading binary data junk is always done with ReadString(), regardless of the session type
% - Added QueryBinaryDataToFile() - use it to transfer files from instrument to PC (see function help for example)
% - Fixed WriteBinaryDataBlock() for SERIAL and SOCKET sessions, added non-mandatory parameter addLFatTheEnd


% Version 1.6.0
% - Added method WriteWithOPC()

% Version 1.5.0
% - Added QueryBinaryDataBlock()
% - Added WriteBinaryDataBlock()

% Version 1.4.0
% - Fixed QueryLongString() method
% - Added variable number of arguments (maximum of 6) as parameters for sprintf to the following methods:
%   Write(), QueryString(), QueryLongString(), QueryDouble(), QueryInteger(), QueryBoolean()
%   Example: instrument.Write('FREQ %0.1f', frequency);

% Version 1.3.0
% - Added settings for SOCKET session type:
%   - termination character enable in constructor
%   - reading binary data junk characters reading is done with Clear()

% Version 1.2.0
% - Added obj.Session.Clear() to the public function obj.ClearStatus()

% Version 1.1.0
% - Error checking changed: Added *STB? at the beginning
% - Added public methods WriteWithSTBpollSync(), QueryWithSTBpollSync(), WriteWithSRQsync(), QueryWithSRQsync()
% - Added private method STBpolling()

% Version 1.0.0
% - first version created
    
    properties
        AddLFtoWriteEnd % Default = false. If true, every string written to the instrument gets linefeed character appended to the end.
        ReadBufferSize % Default = 4096. The buffer size used for reading responses in methods ReadString, QueryString. Methods ReadLongString and QueryLongString use dynamic size that adjust to the required length.
        Session % Interface Ivi.Visa.GlobalResourceManager. Use it to acces the full scope of available interfaces, methods and properties of the Ivi.VISA .NET assembly
        SessionClass % Type of the session (INSTR, SOCKET)
    end
    
    properties (Access=private)
        LastSentCommand % Internal string variable for the purporse of generating exception messages
        VXIcapable % false for SOCKET and SERIAL connections, otherwise true
    end
    
    methods
        function obj = VISA_Instrument(resourceString)
            % Initiates a new VISA connection defined by resourceString
            try 
                assemblyCheck = NET.addAssembly('Ivi.Visa'); %#ok<NASGU>
                %fprintf('Opening Visa .NET assembly...\n');
            catch
                error('Error loading .NET assembly Ivi.Visa');
            end
            
            obj.Session = Ivi.Visa.GlobalResourceManager.Open(resourceString);
            obj.SessionClass = char(obj.Session.ResourceClass);
            
            obj.ReadBufferSize = 4096;
            obj.SetTimeoutMilliseconds(5000);
            
            obj.LastSentCommand = '';
            
            obj.AddLFtoWriteEnd = false;
            
            % VXIcapable is false for socket and serial connections. Otherwise it is true 
            obj.VXIcapable = true;
            if (strcmp(obj.SessionClass, 'SOCKET') == true)
                obj.VXIcapable = false;
            end
            
            if (strcmp(char(obj.Session.HardwareInterfaceType), 'Serial') == true)
                obj.VXIcapable = false;
            end
           
            % for SOCKET and SERIAL sessions:
            % - automatically add LF to every write string
            % - end any reading with receiving LF character
            if (obj.VXIcapable == false)
                obj.AddLFtoWriteEnd = true;
                obj.Session.TerminationCharacterEnabled = 1;
                obj.Session.TerminationCharacter = 10;
            else
                obj.Session.TerminationCharacterEnabled = 0;
                obj.Session.TerminationCharacter = 10;
                obj.Session.Clear();
            end
        end
        
        function Close(obj)
            % Closes the VISA session
            obj.Session.Dispose();
        end
        
        function SetTimeoutMilliseconds (obj, value)
            % Sets timeout for all VISA operations
            % Set it higher than expected operation duration time, e.g. your oscilloscope acquisition time, or your spectrum analyzer sweep time
            %
            % value - timeout value in milliseconds
             obj.Session.TimeoutMilliseconds = value;
        end
        
        function value = GetTimeoutMilliseconds (obj)
            % Returns VISA Timeout.
            %
            % (return) value - current timeout value in milliseconds
            value = obj.Session.TimeoutMilliseconds;
        end
        
        function Write(obj, command, varargin)
            % Sends command string to the instrument
            % If the property AddLFtoWriteEnd is set to true (1), the LF character will be added to the end of the command string
            % Examples:
            % instrument.Write('FREQ:CENT 200E6')
            % instrument.Write('FREQ:CENT %0.1f', frequency)
            % instrument.Write('FREQ%d:CENT %0.1f', window, frequency)
            %
            % command - command string to be sent to the instrument
            % varargin - if the command string contains formatting parameters, list them here. Up to 6 parameters are supported
            if( isempty(command) || (ischar(command)~= 1) )
                throw(MException('VISA_Instrument:Write', 'Parameter "command" is empty or not a string'));
            end
            arguments = length(varargin);
            if (arguments == 0)
                % This makes sure, that even with no arguments, the last '\n' will be replaced by char(10) character (0x0A LF)
                if (length(command) >= 2 && strcmp('\n', command(end-1:end)))
                    command = [ command(1:end-2) char(10) ];
                end
            elseif (arguments == 1)
                command = sprintf(command, varargin{1});
            elseif (arguments == 2)
                command = sprintf(command, varargin{1},varargin{2});
            elseif (arguments == 3)
                command = sprintf(command, varargin{1},varargin{2},varargin{3});
            elseif (arguments == 4)
                command = sprintf(command, varargin{1},varargin{2},varargin{3},varargin{4});
            elseif (arguments == 5)
                command = sprintf(command, varargin{1},varargin{2},varargin{3},varargin{4},varargin{5});
            elseif (arguments == 6)
                command = sprintf(command, varargin{1},varargin{2},varargin{3},varargin{4},varargin{5},varargin{6});
            else
                throw(MException('VISA_Instrument:WriteFormatIntoString','Only up to 6 format parameters are supported'));
            end
                
            obj.LastSentCommand = command;
            
            if (obj.AddLFtoWriteEnd == false)
                obj.Session.RawIO.Write(command)
            else
                obj.Session.RawIO.Write([command char(10)])
            end
        end
        
        % Sends command string to the instrument with added *OPC? query at the end.
        %  Then, it reads the response, and hence waits until the instrument has completed the operation.
        % Examples:
        % instrument.WriteWithOPC('FREQ:CENT 200E6'); % Actual sent string will be 'FREQ:CENT 200E6;*OPC?'
        function WriteWithOPC(obj, command, varargin)
            query = [command, ';*OPC?'];
            obj.QueryString(query, varargin{:});
        end
                
        function response = ReadString(obj)
            % Reads string response from the instrument
            % Potential LF character at the end is trimmed
            % The length of the response must not exceed the current obj.ReadBufferSize
            % If you expect a long string of unknown length, use the ReadLongString() method instead
            % See also READLONGSTRING
            %
            % (return) response - string read from the instrument
            [response, returnStatus] = obj.Session.RawIO.ReadString(obj.ReadBufferSize);
            if (returnStatus == Ivi.Visa.ReadStatus.MaximumCountReached)
                throw(MException('VISA_Instrument:ReadString', ...
                ['Last sent command: "%s"', char(10), ...
                'Not all data from instrument output buffer were read.', char(10), ...
                'The data are probably longer than %d bytes set by the obj.ReadBufferSize property.', char(10), ...
                'Use the QueryLongString or ReadLongString methods.'], obj.LastSentCommand, obj.ReadBufferSize));
            end
            response = char(response);
            response = deblank(response);
        end
        
        function response = ReadLongString(obj)
            % Reads long string response from the instrument
            % Compare to the ReadString() method:
            % - the string length is not limited
            % - the potential LF character at the end is not trimmed.
            %
            % (return) response - string read from the instrument
            % See also READSTRING
            
            response = obj.Session.FormattedIO.ReadString();
            response = char(response);
        end
        
        function response = QueryLongString(obj, query, varargin)
            % Combines the Write() and ReadLongString() into one method
            % Examples:
            % longResponse = instrument.QueryLongString('DISPLAY:WINDOW1:CATALOG?')
            % longResponse = instrument.QueryLongString('DISPLAY:WINDOW%d:CATALOG?', window)
            
            % query - query string to be sent to the instrument
            % varargin - if the query string contains formatting parameters, list them here. Up to 6 parameters are supported
            % (return) response - string read from the instrument
            % See also WRITE, READLONGSTRING
            obj.Write(query, varargin{:});
            response = obj.ReadLongString();
            response = char(response);
        end
        
        function response = QueryString(obj, query, varargin)
            % Combines the Write() and ReadString() into one method
            % Examples:
            % response = instrument.QueryString('FREQ:CENT?')
            % response = instrument.QueryString('FREQ%d:CENT?', window)
            %
            % query - query string to be sent to the instrument
            % varargin - if the query string contains formatting parameters, list them here. Up to 6 parameters are supported
            % (return) response - string read from the instrument. A potential LF character at the end is trimmed
            % See also WRITE, READSTRING
            obj.Write(query, varargin{:});
            response = obj.ReadString();
        end
        
        function double = QueryDouble(obj, query, varargin)
            % Sends a query to the instrument, reads the reponse and converts it to double number
            % Examples:
            % frequency = instrument.QueryDouble('FREQ:CENT?')
            % frequency = instrument.QueryDouble('FREQ%d:CENT?', window)
            %
            % query - query string to be sent to the instrument
            % varargin - if the query string contains formatting parameters, list them here. Up to 6 parameters are supported
            % (return) double - double precision number 
            double = str2double(obj.QueryString(query, varargin{:}));
        end
        
        function integer = QueryInteger(obj, query, varargin)
            % Sends a query to the instrument, reads the reponse and converts it to integer number
            % Examples:
            % frequency = instrument.QueryInteger('FREQ:CENT?')
            % frequency = instrument.QueryInteger('FREQ%d:CENT?', window)
            %
            % query - query string to be sent to the instrument
            % varargin - if the query string contains formatting parameters, list them here. Up to 6 parameters are supported
            % (return) integer - int32 number
            response = obj.QueryString(query, varargin{:});
            integer = str2num(['int32(' response ')']); %#ok<ST2NM>
        end
        
        function boolean = QueryBoolean(obj, query, varargin)
            % Sends a query to the instrument, reads the reponse and converts it to boolean value
            % Examples:
            % output = instrument.QueryBoolean('OUTPUT?')
            % output = instrument.QueryBoolean('OUTPUT%d?', output)
            %
            % query - query string to be sent to the instrument
            % varargin - if the query string contains formatting parameters, list them here. Up to 6 parameters are supported
            % (return) boolean - boolean value
            response = lower(obj.QueryString(query, varargin{:}));
            if strcmp(response, '1') || strcmp(response, 'on') || strcmp(response, 'true')
                boolean = 1;
            else
                boolean = 0;
            end
        end
        
        function listOfDoubles = ReadASCII_ListOfDoubles(obj, maximumNofElements)
            % Reads comma-separated double values from the instrument.
            % The values expected are ACSII-coded. Example '-1.23456E7,+2.98765E+3,+1.1111E+3'
            %
            % maximumNofElements - maximum number of elements to read. If more elements are available, they will not be read out from the instrument
            % (return) listOfDoubles - array of double numbers
            % See also QUERYASCII_LISTOFDOUBLES
            if (~isnumeric(maximumNofElements))
                throw(MException('VISA_Instrument:ReadASCII_ListOfDoubles', 'Variable "maximumNofElements" must be numeric'));
            end
            listOfDoubles = obj.Session.FormattedIO.ReadListOfDouble(maximumNofElements).double;
        end
        
        function listOfDoubles = QueryASCII_ListOfDoubles(obj, query, maximumNofElements, varargin)
            % Combines the Write() and ReadASCII_ListOfDoubles() into one method
            %
            % query - query string to be sent to the instrument
            % maximumNofElements - maximum number of elements to read. If more elements are available, they will not be read out from the instrument
            % varargin - if the query string contains formatting parameters, list them here. Up to 6 parameters are supported
            % (return) listOfDoubles - array of double numbers
            % See also WRITE,READASCII_LISTOFDOUBLES
            obj.Write(query, varargin{:});
            listOfDoubles = obj.ReadASCII_ListOfDoubles(maximumNofElements);
        end
        
        function data = ReadBinaryDataBlock(obj)
            % This method reads data from the instrument that have binary format
            % The binary data start with header - e.g. '#41234' meaning the length of the length is 4 bytes, the length of the binary data is 1234 bytes
            % The following, in our example 1234 bytes are the data content returned by the function
            % Some instruments send an additional linefeed character after the block, this method discards it to prevent further 'Query Interrupted' instrument errors
            %
            % (return) data - array of bytes (uint8) read as binary block from the instrument
            blockLen = obj.ParseBinaryBlockHeader();
            
            if (obj.VXIcapable == false)
                obj.Session.TerminationCharacterEnabled = false;
            end
            
            [data, returnStatus] = obj.Read(blockLen);
            
            if (obj.VXIcapable == false)
                    obj.Session.TerminationCharacterEnabled = true;
            end
                
            % Discard the remaining junk data
            if (returnStatus == Ivi.Visa.ReadStatus.MaximumCountReached)
                obj.ReadString();
            end
        end
        
        function binData = QueryBinaryDataBlock(obj, query, varargin)
            % Combines the Write() and ReadBinaryDataBlock() into one method
            %
            % query - query string to be sent to the instrument
            % varargin - if the command string contains formatting parameters, list them here. Up to 6 parameters are supported
            % (return) binData - binary data returned by the instrument
            obj.Write(query, varargin{:});
            binData = ReadBinaryDataBlock(obj);
        end
        
        function doubleData = ReadBinaryFloatData(obj, swapEndianess)
            % This method reads binary-formatted single-precision data from the instrument and returns them as double-precision array
            % Most of the instruments use single-precision data to return traces or waveforms in a compact 4 bytes/number format
            % This function allows for parsing them and returning an array of double-precision numbers that are commonly used
            % In addition, it allows for swapping Endianess if necessary
            %
            % swapEndianess (non-mandatory, Default = false) If false, no endianess conversion is being done (default behaviour)
            % If the instrument sends data in the Little-endian (most of R&S instruments) and your PC works with Little-endian (most of the PC nowadays), use the swapEndianess=false
            % If the instrument sends data in the Big-endian and your PC works with Little-endian (or vice-versa), use the swapEndianess=true
            % (return) data - array of double numbers read as binary block from the instrument
            data = obj.ReadBinaryDataBlock();
            if nargin < 2
                swapEndianess = false;
            end
            if swapEndianess == false;
                doubleData = typecast(data, 'single');
            else
                data = typecast(data, 'uint32');
                data = swapbytes(data);
                doubleData = typecast(data, 'single');
            end
        end
        
        function doubleData = QueryBinaryFloatData(obj, query, swapEndians)
            % Combines the Write() and ReadBinaryFloatData() into one method
            %
            % query - query string to be sent to the instrument
            % swapEndianess (non-mandatory, Default = false) If false, no endianess conversion is being done (default behaviour)
            % (return) data - array of double numbers read as binary block from the instrument
            % See also WRITE,READBINARYFLOATDATA
            obj.Write(query);
            if nargin < 3
                swapEndians = false;
            end
            doubleData = obj.ReadBinaryFloatData(swapEndians);
        end
        
        function ReadBinaryDataToFile(obj, PCfilePath)
            % This method reads the binary block response from the instrument and saves it to the given file
            %
            % PCfilePath - Path to file in the PC. If the file already exists, it will be overwritten
            fileID = fopen(PCfilePath,'w');
            fwrite(fileID,obj.ReadBinaryDataBlock());
            fclose(fileID);
        end
        
        function QueryBinaryDataToFile(obj, query, PCfilePath)
            % Combines the Write() and ReadBinaryDataToFile() into one method
            % Examples:
            % Copying InstrumentFile.bin from instrument to your PC under the name MyPCfile.bin:
            % QueryBinaryDataToFile('MMEM:DATA? ''InstrumentFile.bin''', 'c:\MyPCfile.bin');
            %
            % query - query string to be sent to the instrument
            % PCfilePath - Path to file in the PC. If the file already exists, it will be overwritten
            obj.Write(query);
            ReadBinaryDataToFile(obj, PCfilePath);
        end
        
        function WriteBinaryDataBlock(obj, command, data, addLFatTheEnd)
            % This method writes a binary data with binary header to the instrument
            % Examples:
            % Copying MyFile.bin file from your PC to the instrument:
            % fileID = fopen('c:\PCfile.bin');
            % data = fread(fileID);
            % fclose(fileID);
            % myInstrument.WriteBinaryDataBlock('MMEM:DATA ''InstrumentFile.bin'',', data);
            %
            % command - command that precedes the binary data block.
            % data - data to be sent to the instrument. The data must be 1xSIZE ot SIZEx1 matrix of uint8 values.
            % addLFatTheEnd - some instrument require LF at the end of the block. In this case, set it to true. SOCKET and SERIAL connections have it automatically set to true. Default value: false
            obj.Session.SendEndEnabled = false;
            
            if nargin < 2
                addLFatTheEnd = false;
            end
            
            if obj.VXIcapable == false
                addLFatTheEnd = true;
            end
               
            % necessary data conversions to uint8
            if isa(data, 'double')
                data = uint8(data);
            end
            
            if ~isa(data, 'uint8')
                throw(MException('VISA_Instrument:WriteBinaryDataBlock', 'Unsupported data type'));
            end
            
            dataSize = size(data);
            size1 = dataSize(1);
            size2 = dataSize(2);
            if size1 == 1
                trueSize = size2;
            elseif size2 == 1
                trueSize = size1;
            else
                throw(MException('VISA_Instrument:WriteBinaryDataBlock', 'Unsupported data size. Only 1 x Size or Size x 1 are allowed.'));
            end
            
            dataSizeString = sprintf('%d',trueSize);
            sizeOfDataSize = length(dataSizeString);
            if sizeOfDataSize < 10
                header = sprintf('%s#%d%s', command, sizeOfDataSize, dataSizeString);
            else
                header = sprintf('%s#(%s)', command, dataSizeString);
            end
            obj.Session.RawIO.Write(header);
            
            if (addLFatTheEnd == false)
                obj.Session.SendEndEnabled = true;
                obj.Session.RawIO.Write(data);
            else
                obj.Session.RawIO.Write(data);
                obj.Session.SendEndEnabled = true;
                obj.Session.RawIO.Write(char(10));
            end
        end
                
        function errors = ReadErrorQueue(obj)
            % Returns cell array of string with all the errors read from the instrument error queue ('*STB?' and 'SYSTEM:ERROR?' query)
            %
            % (return) errors - cell array with all the error in the error queue. If empty, no error was detected
            errors = {};
            stb = obj.QueryInteger('*STB?');
            if bitand(stb, 4) > 0
                while 1
                    response = obj.QueryString('SYST:ERR?');
                    k = strfind(lower(response), '"no error"');
                    if ~isempty (k)
                        break
                    end
                    errors{end+1} = response; %#ok<AGROW>
                end
            end
        end
        
        function ClearStatus(obj)
            % Sends *CLS command and clears the Error Queue if necessary
            obj.Session.Clear();
            obj.QueryString('*CLS;*OPC?');
            obj.ReadErrorQueue();
        end
        
        function ErrorChecking(obj)
            % This method calls the ReadErrorQueue() and if any error is detected, it generates MException()
            %
            % See also READERRORQUEUE
            errors = obj.ReadErrorQueue();
            if ~isempty (errors)
                sizes = size(errors);
                errorsCount = sizes(2);
                allErrors = strjoin(errors, char(10));
                if errorsCount == 1
                    message = 'Instrument reports one error in the error queue';
                else
                    message = sprintf('Instrument reports %d errors in the error queue', errorsCount);
                end
                throw(MException('VISA_Instrument:ErrorChecking', '%s:%s%s', message, char(10), allErrors));
            end
        end
        
        function WriteWithSTBpollSync(obj, command, timeoutms, varargin)
            % This method writes a command with STB polling synchronization mechanism
            %
            % command - command string to be sent to the instrument
            % varargin - if the command string contains formatting parameters, list them here. Up to 6 parameters are supported
            % timeoutms - timeout in milliseconds after which the function throws exception VISA_Instrument:STBpolling
            obj.QueryString('*ESR?');
            fullCmd = sprintf('%s;*OPC', command);
            obj.Write(fullCmd, varargin{:});
            exceptionMessage = sprintf('Writing with OPCsync - Timeout occured. Command: "%s", timeout %d ms', command, timeoutms);
            obj.STBpolling(exceptionMessage, timeoutms);
            obj.QueryString('*ESR?');
        end
        
        function response = QueryWithSTBpollSync(obj, query, timeoutms, varargin)
            % This method sends a query with STB polling synchronization mechanism
            %
            % query - query string to be sent to the instrument
            % timeoutms - timeout in milliseconds after which the function throws exception VISA_Instrument:STBpolling
            % varargin - if the command string contains formatting parameters, list them here. Up to 6 parameters are supported
            % (return) response - string response from instrument
            obj.QueryString('*ESR?');
            fullCmd = sprintf('%s;*OPC', query);
            obj.Write(fullCmd, varargin{:});
            exceptionMessage = sprintf('Querying with OPCsync - Timeout occured. Query: "%s", timeout %d ms', query, timeoutms);
            obj.STBpolling(exceptionMessage, timeoutms);
            response = obj.ReadString();
            obj.QueryString('*ESR?');
        end
        
        function WriteWithSRQsync(obj, command, timeoutms, varargin)
            % This method writes a command with SRQ synchronization mechanism
            %
            % command - command string to be sent to the instrument
            % timeoutms - timeout in milliseconds after which the function throws timeout exception
            % varargin - if the command string contains formatting parameters, list them here. Up to 6 parameters are supported
            obj.QueryString('*ESR?');
            obj.Session.DiscardEvents(Ivi.Visa.EventType.ServiceRequest);
            obj.Session.EnableEvent(Ivi.Visa.EventType.ServiceRequest);
            fullCmd = sprintf('%s;*OPC', command);
            obj.Write(fullCmd, varargin{:});
            obj.Session.WaitOnEvent(Ivi.Visa.EventType.ServiceRequest, timeoutms);
            obj.Session.DisableEvent(Ivi.Visa.EventType.ServiceRequest);
            obj.QueryString('*ESR?');
        end
        
        function response = QueryWithSRQsync(obj, query, timeoutms, varargin)
            % This method sends a query with SRQ synchronization mechanism
            %
            % query - query string to be sent to the instrument
            % timeoutms - timeout in milliseconds after which the function throws timeout exception
            % varargin - if the query string contains formatting parameters, list them here. Up to 6 parameters are supported
            % (return) response - string response from instrument
            obj.QueryString('*ESR?');
            obj.Session.DiscardEvents(Ivi.Visa.EventType.ServiceRequest);
            obj.Session.EnableEvent(Ivi.Visa.EventType.ServiceRequest);
            fullCmd = sprintf('%s;*OPC', query);
            obj.Write(fullCmd, varargin{:});
            obj.Session.WaitOnEvent(Ivi.Visa.EventType.ServiceRequest, timeoutms);
            obj.Session.DisableEvent(Ivi.Visa.EventType.ServiceRequest);
            response = obj.ReadString();
            obj.QueryString('*ESR?');
        end
    end
    
    methods (Access=private)
       
        function [response, returnStatus] = Read(obj, charactersCount)
            % Read method with defined number of characters to be read and the status response
            [response, returnStatus] = obj.Session.RawIO.Read(charactersCount);
            response = response.uint8;
        end
        
        function blockLen = ParseBinaryBlockHeader(obj)
            % Method for parsing the binary block header and returning the parsed length of the binary block
            hashSign = obj.Read(1);
            if char(hashSign) ~= '#'
                nextChars = obj.Read(20);
                throw(MException('VISA_Instrument:Parse_BinaryBlockHeader',...
                'Expecting # (hash) character at the beginning of the binary data block, received "%s%s...."',...
                char(hashSign), char(nextChars)));
            end
            lenOflen = obj.Read(1);
            lenOflen = str2double(char(lenOflen));
            blockLen = obj.Read(lenOflen);
            blockLen = str2double(char(blockLen));
        end
        
        function STBpolling(obj, exceptionMessage, timeoutms)
            % Method for STB polling with progressive polling interval
            timeoutsecs = timeoutms / 1000;
            timerVal = tic;
            while 1
                stb = char(obj.Session.ReadStatusByte());
                                
                k = strfind(stb, 'EventStatusRegister');
                    if ~isempty (k)
                        break
                    end
                
                elapsed = toc(timerVal);
                if elapsed > timeoutsecs
                    throw(MException('VISA_Instrument:STBpolling', exceptionMessage));
                else
                    % Progressive delay
                    if elapsed > 0.01
                        if elapsed > 0.1
                            if elapsed > 1
                                if elapsed > 10
                                    pause(0.5)
                                else
                                    pause(0.1)
                                end
                            else
                                pause(0.01)
                            end
                        else
                            pause(0.001)
                        end
                    end
                end
            end     
        end
    end
end


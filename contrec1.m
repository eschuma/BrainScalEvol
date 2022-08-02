%
%   contrec1.m
%       Continuous 1-channel A/D recording
%   
%   For use with Tucker-Davis Technologies System 3, using an RX8
%   multiprocessor
%
%   CONNECTIONS:
%       ADC1 ... from recording amplifier
%

%% Get recording parameters
disp(' ')
ampgain = input('Amplifier Gain? ');
disp(' ')
if (isempty(ampgain))
    error('Enter Amplifier Gain!!!')
end
recdur = input('Recording Duration? (s) ');
if (isempty(recdur))
    error('Enter Recording Duration!!!')
end
numreps = input('Number of Reps? ');
if (isempty(recdur))
    error('Enter Number of Reps!!!')
end
disp(' ')

%% Load circuit and start to run
RX8 = actxcontrol('RPco.x',[5 5 26 26]);
RX8.ConnectRX8('GB', 1);
if isdir('N:\')==1,
    RX8.LoadCOF('N:\Matlab\contrec1.rcx');
else
    RX8.LoadCOF('C:\Documents and Settings\Carlson Lab\My Documents\Matlab\contrec1.rcx');
end
srate = RX8.GetSFreq;
RX8.Run;
pause(1)    % pausing to avoid effect of transient noise at circuit startup
maxdur = 1e6/srate;
if recdur>maxdur,
    error('Recording Duration is Too Long!!!')
end

%% Set recording parameters
npts = ceil(recdur*srate);                 % number of points in sample
rectime = (1/srate)*[1:1:npts];             % recording time points

%% Check for errors
status = double(RX8.GetStatus);     % Gets the status
if bitget(status,1)==0;             % Checks for connection
    disp('Error connecting to RX8'); return;
elseif bitget(status,2)==0;     % Checks for errors in loading circuit
    disp('Error loading circuit'); return;
elseif bitget(status,3)==0      % Checks for errors in running circuit
    disp('Error running circuit'); return;
end

%% Record Data
if all(bitget(RX8.GetStatus,1:3))
    
    for i=1:numreps,
    
        % Send buffer information to RX8
        RX8.SetTagVal('schmittsize',npts);
        RX8.SetTagVal('bufsize',npts+1);

        % Trigger recording
        RX8.SoftTrg(1);

        % Wait until buffer is full
        curindex=RX8.GetTagVal('curindex');
        while(curindex<npts)
            curindex=RX8.GetTagVal('curindex');
        end
        bufdata=RX8.ReadTagV('bufdata', 0, npts);
        if i==1,
            recdata = bufdata./ampgain;
        else
            recdata = [recdata; bufdata./ampgain];
        end
        
    end

    % Stop
    RX8.Halt;
    
end

%% Display Data
figrec = figure('Name','Continuous Recording','NumberTitle','off','Position',[50 150 1200 600],'Color','w');
plot(rectime,recdata)
xlim([min(rectime) max(rectime)])
xlabel('Time (s)')
ylabel('Potential (V)')
zoom on

%% Save Data
disp(' ')
savdata = input('Save Data to File? (y/n) ','s');
    if savdata=='y',
    [filename,pathname] = uiputfile('*.crc1.mat','Enter file name');
    if pathname == 0
        beep
        return
    end
    save (fullfile(pathname,filename), 'srate', 'rectime', 'recdata');
    cd(char(pathname))
end
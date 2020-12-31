clc;clear

% Initialise the library
CEDS64ML_path = '../CEDS64ML/'; % the location of CED64ML library on this computer
addpath(CEDS64ML_path); % add CED64ML to the Matlab search path
setenv('CEDS64ML', CEDS64ML_path); % set environmental variable CED64ML
CEDS64LoadLib( CEDS64ML_path ); % load ceds64int.dll

% load PurkinjeSort_challange_dataset
load('PurkinjeSort_challange_dataset.mat')

% open smrx file
fhand = CEDS64Create( './PurkinjeSort_challange_dataset.smr' );
if (fhand <= 0); CEDS64ErrorMessage(fhand); unloadlibrary ceds64int; return; end

filesize1 = CEDS64FileSize(fhand);
maxtime1 = CEDS64MaxTime(fhand);

% set and get file time base
CEDS64TimeBase( fhand, 1/sample_rate );
timebase = CEDS64TimeBase( fhand );

version = CEDS64Version(fhand);

wavevecdoubl = ch_data;

% create real wave channel
wavechan2 = CEDS64GetFreeChan( fhand );
createret = CEDS64SetWaveChan( fhand, wavechan2, 1, 9, sample_rate );
if createret ~= 0, warning('realwave channel 1 not created correctly'); end;
CEDS64ChanTitle( fhand, wavechan2, 'Data');
CEDS64ChanComment( fhand, wavechan2, 'PurkinjeSort (Psort) chalenge');
CEDS64ChanUnits( fhand, wavechan2, 'uV' );

% fill the RealWave channel
sTime = CEDS64SecsToTicks( fhand, 0 ); % offset start by 0 seconds
fillret = CEDS64WriteWave( fhand, wavechan2, wavevecdoubl, sTime );
if fillret < 0, warning('RealWave channel 1 not filled correctly'); end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% close the file
CEDS64CloseAll();
% unload ceds64int.dll
unloadlibrary ceds64int;




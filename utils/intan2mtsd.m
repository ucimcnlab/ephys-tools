function [eeg,par] = intan2mtsd(fnbin,varargin)
%
% Reads raw data from the .bin file and converts to mV.  Data is returned
% as Mtsd object.  Raw data file and settings file are required (.bin and
% .xml)
%
% required input:
% fnbin - file name of binary file
%
% optional name,value args:
% 'channels' - which channels to plot, default: [] (empty=all)
% 'epoch' - time range to plot in sec: [start end], default: [] (all)
% 'doread' - true/false if you really want to read the file.  Otherwise,
%            information about the file is returned. default: true
% 'settingsfile' - name of xml settings file, default: 'settings.xml'
%
% Output
% eeg - data as a Mtsd object
% par - structure with file information and parameters. Fields are:
%          binaryfile - name of binary file
%       setttingsfile - name of settings file
%       validchannels - all valid channels on the probe
%            channels - requested channels to read
%               epoch - time interval specified for reading
%                  fs - sampling frequency of the data
%               scale - multiplier used to convert data to volts
%            nsamprec - number of data samples in the file
%              recdur - duration of recording in sec
%              doread - true/false indicating if read was requested
%

par = inputParser;
addParameter(par, 'channels', [], @isnumeric);
addParameter(par, 'epoch',[]);
addParameter(par, 'doread',true , @islogical);
addParameter(par, 'settingsfile', 'settings.xml', @(x)ischar(x)||isstring(x));
parse(par, varargin{:});
par = par.Results;

% validate inputs and parameters
if ~isfile(fnbin)
    error('file not found: %s',fnbin)
end
par.binaryfile = fnbin;
if ~isfile(par.settingsfile)
    error('file not found: %s',par.settingsfile)
end

hdr = intan_settings(par.settingsfile);
fs = hdr.fs;
scale = hdr.AnalogScaleVolts;
par.fs = fs;
par.scale = scale;
if contains(fnbin,'hc','ignorecase',true)
    ch = hdr.nchan(3:4);
else
    ch = hdr.nchan(1:2);
end
validchans = [ch(1).chans(ch(1).okchan) ch(2).chans(ch(2).okchan)+128]+1;
nvalid = length(validchans);
par.validchans = validchans;
if ~isempty(par.channels)
    chans = par.channels;
    ix = ismember(chans,validchans);
    if ~all(ix)
        s = num2str(chans(~ix));
        error('invalid channels: %s',s)
    end
else
    chans = validchans;
end
f = dir(fnbin);
nsamprec = f.bytes/nvalid/2;
recdur = nsamprec/fs;
par.nsamprec = nsamprec;
par.recdur = recdur;
if ~isempty(par.epoch)
    epoch = par.epoch;
    if length(epoch)~=2 || epoch(1)<0 || (~isinf(epoch(2))&&epoch(2)>recdur)
        error('invalid epoch specified')
    end
    if isinf(epoch(2))
        epoch(2) = recdur;
    end
else
    epoch = [0 recdur];
end

% read in data
if par.doread
    ix = round(epoch*fs);
    nsamp = diff(ix)+1;
    ncr = length(chans);
    mem = memory;
    mx = mem.MaxPossibleArrayBytes;
    nbfile = ncr*nsamp*2;
    nbmat = ncr*nsamp*4;
    if nbfile+nbmat > mx
        warning('not enough memory')
        r = input('continue with read??? (y/[n]) ','s');
        if ~strcmpi(r,'y')
            return
        end
    end
    offset = round(ix(1)*nvalid*2) + (chans(1)-1)*2;
    fid = fopen(fnbin,'rb');
    fseek(fid,offset,'bof');
    a = fread(fid,[nvalid,nsamp],'int16=>int16')';
    fclose(fid);
    nsamp = size(a,1);
    a = single(a(:,chans)) .* scale;
    dt = 1/fs;
    t = 0:dt:dt*(nsamp-1);
    t = t + epoch(1);
    eeg = Mtsd(t,a);
end
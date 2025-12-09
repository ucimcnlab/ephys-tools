function eeg = npy2mtsd(fname,varargin)
%
% Reads a numpy LFP and returns a Mtsd object
%
% fname - numpy file name containing LFP data
%
% Optional name/value pairs
% 'fs' - sampling frequency (Hz). If not given, tries to find timestamp file.
% 'channel' - optional single channel to extract, default is all
% 'chunk' - optional start/stop indices for import, default is whole file
%

par = inputParser;
addParameter(par, 'fs', [], @(x)isnumeric(x)&&length(x)<2)
addParameter(par, 'channel', [], @(x)isnumeric(x)&&length(x)<2)
addParameter(par, 'chunk', [], @(x)isnumeric(x)&&length(x)==2)
parse(par, varargin{:});
par = par.Results;


shape = readNPYheader(fname);
if isscalar(shape)
    shape = [shape 1];
end
if isempty(par.chunk)
    ix = [1 max(shape)];
else
    ix = par.chunk;
end
ic = par.channel;
fs = par.fs;

if isempty(fs)
    i = strfind(fname,'.npy')-1;
    fn = [fname(1:i) 'ts.npy'];
    if isfile(fn)
        x = readnpychunk(fn,[1 2]);
        fs = 1 / (x(2)-x(1));
        if round(fs) == fs
            fmt = '%d';
        else
            fmt = '%.2g';
        end
        fprintf('\nFound timestamp file: %s\n',fn);
        fprintf(['Calculated sampling frequency: ' fmt ' Hz\n'],fs);
    else
        error('need sampling frequency or timestamp file: %s',fn)
    end
end

if ~isempty(ic) && isscalar(ic)
    a = readnpychan(fname,ic);
else
    a = readnpychunk(fname,ix);
end
if size(a,2) > size(a,1)
    a = a';
end
[nsamp,nch] = size(a);
dt = 1/fs;
t = (0:dt:dt*(nsamp-1))' + ix(1)*dt;
eeg = Mtsd(t,a);
1;
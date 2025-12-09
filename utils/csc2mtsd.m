function eeg = csc2mtsd(fname,ts)
%
% Read neuralynx csc file (LFP data) and return Mtsd object.  
%
% fname - name of csc file
% ts - optional start/end timestamps (seconds) to extract. 
%      default: [] (whole file)
%

% setup nlx input args
fs = [1 0 1 1 1]; % get timestamps, sampling frequency, # valid samples, data
hd = 1; % get header
if nargin < 2 || isempty(ts)
    exmode = 1;
    exmodev = [];
else
    exmode = 4;
    exmodev = round(ts*1e6); % convert sec to usec
end

% read in data and get parameters
[t,~,valid,x,hdr] = Nlx2MatCSC(fname,fs,hd,exmode,exmodev);
ix = contains(hdr,'ADBitVolts');
line = strsplit(hdr{ix});
scale = str2double(line{2})*1000;
dt = (t(2)-t(1)) / 512 / 1e6; % time between samples in sec
if ~any(valid)
    valid(1:end) = 512;
end

% construct signal by checking num valid samples in each record, i.e don't
% assume that signal is continuous.
nrec = length(t);
len = nrec*512;
tt = zeros(1,len);
ss = zeros(1,len);
pos = 1;
for i = 1:nrec
    tcur = double(t(i))/1e6;
    nv = valid(i);
    ix = pos:pos+nv-1;
    ss(ix) = x(1:nv,i)';
    tt(ix) = tcur:dt:tcur+(nv-1)*dt;
    pos = pos + nv;
end
pos = pos-1;
eeg = Mtsd(tt(1:pos)',ss(1:pos)'*scale);
function v = vt2mtsd(varargin)
%
% Read in neuralynx tracker file and make Mtsd object.  Extracts the 
% timestamps, x and y coordinates from the file.  The 'voltage' field
% of the Mtsd is a npointx2 matrix of x,y coordinates.  By default,
% missing data points (x==0 & y==0) are linearly interpolated.
%
% Optional name/value pair arguments
% 'fname' - optional name of tracker file, default: VT1.nvt
% 'epoch' - optional 1x2 vector of start/end timestamps for extracting.
%           default is [] and the whole file is read.
% 'interp' - true/false to interpolate missing data. default: true
%
% Returns
% v - Mtsd of timestamps (seconds), and xy coordinates
%


par = inputParser;
addParameter(par, 'fname', 'VT1.nvt', @isfile);
addParameter(par, 'epoch', [], @(x)numel(x)==2&&isnumeric(x));
addParameter(par, 'interp', true, @islogical);
parse(par, varargin{:});

fname = par.Results.fname;
tr = par.Results.epoch;
dointerp = par.Results.interp;

if isempty(tr)
    exmode = 1; % extract all
else
    exmode = 4; % extract timestamp range
end

gethdr = 0;
% extract timestamps, x, y.
fldsel = [1 1 1 0 0 0];

[t,x,y] = Nlx2MatVT(fname,fldsel,gethdr,exmode,tr);
t = t/1e6;
xy = [x;y]';

if dointerp
    ix = xy(:,1)==0 & xy(:,2)==0;
    tbad = t(ix);
    tgood = t(~ix);
    xygood = xy(~ix,:);
    xyi = interp1(tgood,xygood,tbad);
    xy(ix,:) = xyi;
end

v = Mtsd(t,xy);

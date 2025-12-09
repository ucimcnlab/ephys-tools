function maps = ratemap1d(tbin,vel,spikes,varargin)
%
% Calculate the 1D ratemap
%
% Required Inputs:
% tbin - ntrial x nbin+1 matrix of time bin edges for binning spikes
% vel - ntrial x nbin matrix of running speeds
% spikes - either a cell array of spiketimes or a S matrix of spike times
%
% Optional name,value arguments:
% 'velocitythresh' - ignore bins where speed is below threshold. default: 2
% 'sdsmooth' - standard deviation in bins for gaussian smoothing. default: 1
%
% Output:
% maps structure with fields:
%  ocmap - occupancy map in seconds (ntrial x nbin)
%  ocmapsm - smoothed occupancy
%  ocmap1d - trial averaged occupancy (1 x nbin)
%  ocmap1dsm - smoothed average occupancy
%
%  spkmap - histogram of spike counts (ntrial x nbin)
%  spkmapsm - smoothed spike map
%  spkmap1d - trial averaged spike map (1 x nbin)
%  spkmap1dsm - smoothed average spike map
%
%  rmap - rate map (spike counts normalized by occupancy, ntrial x nbin)
%  rmapsm - smoothed rate map
%  rmap1d - trial averaged rate map (1 x nbin)
%  rmap1dsm - smoothed average rate map
%

par = inputParser;
addParameter(par, 'velocitythresh', 2, @(x)isscalar(x)&&isnumeric(x));
addParameter(par, 'sdsmooth', 1, @(x)isscalar(x)&&isnumeric(x));
parse(par, varargin{:});
par = par.Results;

velthresh = par.velocitythresh;
win = par.sdsmooth*5; % 'smoothdata' sd is 1/5 of specified window size

if isa(spikes,'cell')
    ncell = length(spikes);
elseif isa(spikes,'Msmat')
    ncell = spikes.ncell;
    spikes = spikes.cells;
else
    error('spikes argument must be a cell array or a Msmat object')
end

oc = diff(tbin,1,2);
if any(oc<0,'all')
    error('found negative occupancy, need to correct time bins')
end
[nlap,nbin] = size(oc);

% count spikes to make spike map
ed = tbin';
ed = ed(:);
c = cellfun(@(x)histcounts(x.t,ed),spikes,'uniformoutput',false);
spkmap = vertcat(c{:})';
spkmap(nbin+1:nbin+1:end,:) = [];
spkmap = reshape(spkmap,[nbin nlap ncell]);
spkmap = permute(spkmap,[2 1 3]);

% ignore bins if mouse is not running
ix = vel<velthresh;
oc(ix) = 0;
ix = repmat(ix,[1 1 ncell]);
spkmap(ix) = 0;

% make various versions of ratemap
rmap = bsxfun(@rdivide,spkmap,oc);
rmap(isnan(rmap)|isinf(rmap)) = 0;
ocsm = smoothdata(oc,2,'gaussian',win);
spkmapsm = smoothdata(spkmap,2,'gaussian',win);
rmapsm = bsxfun(@rdivide,spkmapsm,ocsm);
rmapsm(isnan(rmapsm)|isinf(rmapsm)) = 0;

% make output struct
maps.ocmap = oc;
maps.ocmapsm = ocsm;
maps.ocmap1d = mean(oc,'omitnan');
maps.ocmap1dsm = smoothdata(maps.ocmap1d,'gaussian',win);

maps.spkmap = spkmap;
maps.spkmapsm = spkmapsm;
maps.spkmap1d = squeeze(mean(spkmap,'omitnan'))';
maps.spkmap1dsm = smoothdata(maps.spkmap1d,2,'gaussian',win);

maps.rmap = rmap;
maps.rmapsm = rmapsm;
maps.rmap1d = bsxfun(@rdivide,maps.spkmap1d,maps.ocmap1d);
maps.rmap1dsm = bsxfun(@rdivide,maps.spkmap1dsm,maps.ocmap1dsm);
maps.bins = tbin;
1;
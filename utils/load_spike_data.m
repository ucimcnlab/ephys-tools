function dat = load_spike_data(varargin)
%
% Loads spike data from Kilosort and optionally classifies clusters as
% pyramidal or interneuron based on firing rate and waveform
%
% optional name,value arguments:
% 'fn_info' - name of kilosort cluster info file. default: 'cluster_info.tsv'
% 'fn_clu' - name of kilosort spike cluster file. default: 'spike_clusters.npy'
% 'fn_times' - name of kilosort spike times file. default: 'spike_times.npy'
% 'fs' - sampling frequency in Hz. default: 30,000
% 'groups' - which group labels to extract.  Can be a single group or a 
%            cell array. default: 'good'
% 'epoch' - limit data extraction to time interval as a [1 2] vector of
%           start/end times in seconds.  Default is to use global min/max
%           of all spiketimes pooled across cells.
% 'doclassify' - true/false if you want to classify clusters as pyramidal
%                or interneuron based on peak-trough duration and firing
%                rate.  Requires bombcell waveform duration file.
% 'class_fr' - firing rate threshold for classification in Hz. default: 8
% 'class_dur' - waveform peak to trough duration for classification in
%               microseconds. default: 400
% 'fn_wav' - name of bombcell waveform duration file.
%            default: 'cluster_waveform_duration.tsv'
%
% Outputs structure 'dat' with fields:
%         clu: cluster number
%       depth: position of cell along electrode (microns)
%          ch: channel of max amplitude
%  spiketimes: spiketimes in seconds
%          fr: average firing rate (Hz)
%    celltype: 'pyramidal' or 'interneuron' (if requested)
%

par = inputParser;
addParameter(par, 'fn_info', 'cluster_info.tsv', @(x)ischar(x)||isstring(x));
addParameter(par, 'fn_clu', 'spike_clusters.npy', @(x)ischar(x)||isstring(x));
addParameter(par, 'fn_times', 'spike_times.npy', @(x)ischar(x)||isstring(x));
addParameter(par, 'fn_wav', 'cluster_waveform_duration.tsv', @(x)ischar(x)||isstring(x));
addParameter(par, 'fs', 30000, @(x)isscalar(x)&&isnumeric(x));
addParameter(par, 'groups', 'good', @(x)ischar(x)||iscell(x));
addParameter(par, 'epoch', [], @(x)isnumeric(x)&&numel(x)==2);
addParameter(par, 'doclassify', false, @islogical);
addParameter(par, 'class_fr', 8, @(x)isscalar(x)&&isnumeric(x));
addParameter(par, 'class_dur', 40, @(x)isscalar(x)&&isnumeric(x));
parse(par, varargin{:});
fninfo = par.Results.fn_info;
fnclu = par.Results.fn_clu;
fntimes = par.Results.fn_times;
fnwav = par.Results.fn_wav;
doclassify = par.Results.doclassify;
frclass = par.Results.class_fr;
durclass = par.Results.class_dur;
fs = par.Results.fs;
grps = par.Results.groups;
epoch = par.Results.epoch;


T = readtable(fninfo,'filetype','text','delimiter','\t');
cluall = readNPY(fnclu);
tall = readNPY(fntimes);
tall = double(tall) ./ fs;

if iscell(grps)
    nc = size(T,1);
    ix = false(nc,1);
    for i = 1:length(grps)
        ix = ix | strcmpi(T.group,grps{i});
    end
else
    ix = strcmpi(T.group,grps);
end

ngood = nnz(ix);
dat.clu = T.cluster_id(ix);
dat.depth = T.depth(ix);
dat.ch = T.ch(ix)+1;
dat.spiketimes = cell(1,ngood);
for i = 1:ngood
    ix = cluall==dat.clu(i);
    dat.spiketimes{i} = tall(ix);
end
s = Msmat(dat.spiketimes);
if ~isempty(epoch)
    s = s.restrict(epoch);
end
dat.fr = s.nspike' ./ diff(s.tlim);
dat.spiketimes = s.cells;

if doclassify
    D = readtable(fnwav,'filetype','text','delimiter','\t');
    dat.celltype = cell(1,ngood);
    for i = 1:ngood
        ix = D.cluster_id == dat.clu(i);
        if nnz(ix)==1
            dur = D.waveform_duration(ix);
            if dur > 400
                if dat.fr(i) < 8
                    dat.celltype{i} = 'pyramidal';
                else
                    dat.celltype{i} = 'interneuron';
                end
            else
                dat.celltype{i} = 'interneuron';
            end
        end
    end
end
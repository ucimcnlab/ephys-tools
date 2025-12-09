function axall = ratemap_figs(rmap,varargin)
%
% Make ratemap images for a number of cells.  The data 'rmap' should be 
% organized ntrial x nbin x ncell.  Makes a 1xn row of axes in 
% the figure window, where n can be specified.  Generates as many figures
% as necessary depending on number of cells.  Optionally plots the trial
% averaged ratemap below each trial ratemap.
%
% required input:
% rmap - the ratemaps to plot.  Should be ntrial x nbin x ncell
%
% optional name/value pairs:
% 'cells_per_figure' - number of plots per figure. default: 6
% 'fig_width' - figure width. default: 1200
% 'fig_height' - figure height. default: 400
% 'do_avg' - true/false to plot average. default: true
% 'cmap' - colormap to use. default: a monochrome blue-white map
% 'clusterid' - vector of numeric cluster ID's to be used for figure titles.
%               If not given, figures are numbered from 1:ncell.
%
% Returns a cell array [1,nfig] of axes handles
%

par = inputParser;
addParameter(par, 'cells_per_figure', 6, @(x)isscalar(x)&&isnumeric(x));
addParameter(par, 'fig_width', 1200, @(x)isscalar(x)&&isnumeric(x));
addParameter(par, 'fig_height', 400, @(x)isscalar(x)&&isnumeric(x));
addParameter(par, 'do_avg', true, @islogical);
addParameter(par, 'cmap', [], @(x)isnumeric(x)&&isequal(shape(x),[256 3]));
addParameter(par, 'clusterid', [], @isnumeric);

parse(par, varargin{:});

cpf = par.Results.cells_per_figure;
fw = par.Results.fig_width;
fh = par.Results.fig_height;
doavg = par.Results.do_avg;
cmap = par.Results.cmap;
if isempty(cmap)
    cmap = [linspace(.93,0,256); linspace(.96,.33,256); linspace(1,.65,256)]';
end
clus = par.Results.clusterid;

[ntrial,nbin,ncell] = size(rmap);
if isempty(clus)
    clus = 1:ncell;
end
if length(clus) ~= ncell
    error('number of cluster IDs and cells does not match')
end
nfig = ceil(ncell/cpf);
axall = cell(1,nfig);
ic = 0;
for ifig = 1:nfig
    figure('position',[400 400 fw fh]);
    o = tiledlayout(10,cpf,'padding','compact','tilespacing','compact');
    ax = gobjects(2,cpf);
    if ifig < nfig
        nc = cpf;
    else
        nc = mod(ncell,cpf);
    end
    for i = 1:nc
        ic = ic + 1;
        ax(1,i) = nexttile([9 1]);
        a = rmap(:,:,ic);
        imagesc(a)
        mx = prctile(a(:),99.5);
        clim([0 mx])
        av = mean(a);
        fr = mean(av);
        st = sprintf('clu: %d,  %.1f Hz',clus(ic),fr);
        title(st)

        if doavg
            ax(2,i) = nexttile(9*cpf+i);
            imagesc(av)
            clim([0 mx])
            if i==1
                ylabel('mean')
            end
        end
    end
    set(ax(1,1:nc),'xtick',[])
    if doavg
        set(ax(2,1:nc),'ytick',[])
    end
    o.XLabel.String = 'Bins';
    o.YLabel.String = 'Trials';
    axall{ifig} = ax;
    colormap(cmap)
end
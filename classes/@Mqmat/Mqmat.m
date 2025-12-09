classdef Mqmat
    % Q-matrix class representing matrix of spike counts based on discrete
    % time bins.  Matrix is  nbin x ncell and must be continuous in time
    % (unlike the Mtsd class, which can have gaps)

    properties (GetAccess = public, SetAccess = private)
        tstart (1,1) double % start time of q-matrix (seconds)
        dt (1,1) double % bin size in seconds
        data double % binned spike counts (nbin x ncell)
        tlim (1,2) double % start/end times for q-matrix (seconds)
        ncell (1,1) double % number of cells in q-matrix
    end


    methods (Access = public)
        function q = Mqmat(varargin)
            % constructor - 2 ways to make q matrix
            % 1. From raw data.  Input arguments must be:
            %    tstart - start time (seconds)
            %    dt - time difference (seconds) between bins
            %    data - matrix of binned data, shape should be nbin x ncell,
            %           user is warned if ncell > nbin since this is unusual
            %
            % 2. From s-matrix (cell array of timestamps)
            %    Input arguments are:
            %    s - s matrix of timestamps in seconds
            %    binsize - binsize (dt) in seconds
            %    tlim - Optional 2-element vector of start/end timestamps
            %           in seconds.  If not given, start/end times are
            %           determined by min/max spike times
            switch class(varargin{1})
                case 'double'
                    if nargin < 3
                        error('arguments should be: t0, dt, data')
                    end
                    q.tstart = varargin{1};
                    q.dt = varargin{2};
                    q.data = varargin{3};
                    q.tlim = q.calctlim;
                    [nr,nc] = size(q.data);
                    q.ncell = nc;
                    if nc > nr
                        warning('Number of cells is greater than number of time bins!!')
                    end

                case {'cell','Msmat'}
                    if nargin < 2
                        error('need bin size')
                    end
                    s = varargin{1};
                    if isa(s,'cell')
                        s = Msmat(s);
                    end
                    dt = varargin{2};
                    if nargin==3 && length(varargin{3})==2
                        tlim = varargin{3};
                        s = s.restrict(tlim);
                    else
                        tlim = s.tlim;
                    end

                    % build row,column indices for the location of
                    % spikes in the q-matrix
                    nspiketot = sum(s.nspike);
                    icell = zeros(nspiketot,1);
                    itime = zeros(nspiketot,1);
                    icur = 1;
                    for i = 1:s.ncell
                        if s.nspike(i) > 0
                            ix = icur:icur+s.nspike(i)-1;
                            icell(ix) = zeros(s.nspike(i),1) + i;
                            itime(ix) = s.cells{i}.t;
                            icur = icur + s.nspike(i);
                        end
                    end
                    % convert spike times to bin indices
                    itime = floor((itime - tlim(1))/dt) + 1;
                    nbin = floor(diff(tlim)/dt);
                    % discard any spikes past end of q-matrix
                    ix = itime>nbin;
                    itime(ix) = [];
                    icell(ix) = [];
                    nspiketot = nspiketot - nnz(ix);

                    % generate sparse nbin x ncell matrix of spike counts
                    if isempty(itime)
                        dat = sparse(nbin,s.ncell);
                    else
                        n = ones(nspiketot,1);
                        dat = sparse(itime,icell,n,nbin,s.ncell);
                    end
                    q.tstart = tlim(1);
                    q.dt = dt;
                    q.data = dat;
                    q.tlim = q.calctlim;
                    q.ncell = s.ncell;

                otherwise
                    error('input must be Smatrix, cell array of spike times, or raw data as (tstart,dt,data)')
            end
        end

        function ed = binedges(q)
            % get the times of the bin edges for the q-matrix
            % ed = binedges(q)
            %
            % ed - vector of bin edge times (seconds)
            n = size(q.data,1);
            ed = (0:q.dt:q.dt*n) + q.tstart;
        end

        function cn = bincenters(q)
            % get the times of the bin centers for the q-matrix
            % cen = bincenters(q)
            % 
            % cen - vector of bin center times (seconds)
            a = q.tstart + q.dt/2;
            n = size(q.data,1)-1;
            cn = (0:q.dt:q.dt*n) + a;
        end

        function r = restrict(q,tr)
            % restricts q-matrix to start/end times in 'tr' 
            % r = restrict(q,tr)
            % 
            % tr - n x 2 start/end times (seconds)
            % 
            % r - restricted q-matrix
            if numel(tr) ~= 2
                error('restrict time must contain 2 values')
            end
            ix = floor((tr - q.tstart)/q.dt) + 1;
            tstartnew = (ix(1)-1)*q.dt + q.tstart;
            r = Mqmat(tstartnew, q.dt, q.data(ix(1):ix(2),:));
        end

        function image(q,tp)
            % makes image of the q-matrix with optional timestamp patches
            % image(q,tp)
            %
            % tp - optional n x 2 start/end times (seconds) for patches
            x = q.bincenters;
            y = 1:q.ncell;
            figure;
            imagesc(x,y,q.data')
            xlabel('Time (sec)')
            ylabel('Cells')
            set(gca,'position',[.05 .11 .93 .82])
            if nargin==2 && ~isempty(tp)
                tp = restrict_ts(tp,q.tlim);
                x = [tp fliplr(tp)]';
                y = [ylim;ylim];
                y = repmat(y(:),1,size(x,2));
                patch(x,y,'w','edgecolor','w','facealpha',.35,'linewidth',1)
            end
        end

        function [av,t,se] = trigavg(o,tr,dur,noavg)
            arguments (Input)
                o (1,1) Mqmat
                tr (:,1) double
                dur (1,1) double
                noavg (1,1) logical = false
            end
            ix = tr > o.tlim(1)+dur & tr < o.tlim(2)-dur;
            tr = tr(ix);
            n = numel(tr);
            nbin = round(dur/o.dt);
            ed = o.binedges;
            a = zeros(o.ncell,2*nbin,n);
            for i = 1:n
                ix = binsearch(ed,tr(i));
                a(:,:,i) = o.data(ix-nbin:ix+nbin-1,:)';
            end
            if noavg
                av = a;
            else
                av = mean(a,3);
            end
            if nargout>1
                t = (-nbin:nbin-1)*o.dt+o.dt/2;
                if nargout>2
                    se = std(a,[],3) ./ sqrt(n);
                end
            end
        end
    end


    methods (Access = private)
        function tlim = calctlim(q)
            e = q.tstart + size(q.data,1)*q.dt;
            tlim = [q.tstart e];
        end
    end
end

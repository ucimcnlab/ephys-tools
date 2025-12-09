classdef Msmat
% class for a cell array of Mts objects, called a S-matrix.  
% A Mts is a vector of spike times for a single cell, 
% Msmat is the spike times for a number of cells

    properties (GetAccess = public, SetAccess = private)
        cells (1,:) cell % cell array holding individual Mts objects
        ncell (1,1) double % number of cells
        tlim (1,2) double % global time limits of all cells
        nspike (1,:) double % vector of spike counts
    end

    
    methods (Access = public)
        function s = Msmat(c)
        % constructor. Makes S-matrix from spike times.
        % s = Msmat(c)
        % 
        % c -  cell array of spike times (in seconds) or Mts objects 
        % 
        % s - S matrix
            arguments (Input)
                c (1,:) cell
            end

            s.ncell = length(c);
            for i = 1:s.ncell
                if isa(c{i},'Mts')
                    t = c{i}.t;
                else
                    t = c{i};
                end
                c{i} = Mts(t(:));
            end
            s.cells = c;
            s.nspike = s.calcnspike;
            s.tlim = s.calctlim;
        end

         function s = restrict(s,tr)
         % restricts s-matrix to times in 'tr'.
         % s = restrict(s,tr)
         % 
         % tr - n x 2 start/end timestamps (seconds)
         %
         % s - S matrix restricted to times in 'tr'
            arguments (Input)
                s (1,1) Msmat
                tr (:,2) double
            end

            for i = 1:s.ncell
                c = s.cells{i};
                s.cells{i} = c.restrict(tr);
            end
            s.nspike = s.calcnspike;
            s.tlim = s.calctlim;
         end

         function raster(s,tp)
         % raster plot of spikes, with optional timestamp patches.  
         % raster(s,tp)
         % 
         % tp - optional nx2 start/end times (seconds) to patch
             figure; 
             axes('position',[.05 .11 .9 .82])
             hold on
             for i = 1:s.ncell
                 plot(s.cells{i}.t,zeros(1,s.nspike(i))+i,'k|','markersize',16,'linewidth',1.5)
             end
             ylim([0 s.ncell+1])
             xlabel('Time (sec)')

             if nargin==2 && ~isempty(tp)
                tp = restrict_ts(tp,s.tlim);
                x = [tp fliplr(tp)]';
                y = [ylim;ylim];
                y = repmat(y(:),1,size(x,2));
                patch(x,y,'r','edgecolor','none','facealpha',.25);
            end
         end
    end

    methods (Access = private)
        function ns = calcnspike(s)
        % number of spikes for each cell
            ns = zeros(1,s.ncell);
            for i = 1:s.ncell
                ns(i) = numel(s.cells{i}.t);
            end
        end

        function tlim = calctlim(s)
        % min/max spike times of entire s-matrix
            mx = -inf;
            mn = inf;
            for i = 1:s.ncell
                if ~isempty(s.cells{i}.t)
                    mx = max(mx, max(s.cells{i}.t));
                    mn = min(mn, min(s.cells{i}.t));
                end
            end
            tlim = [mn mx];
        end
    end
end
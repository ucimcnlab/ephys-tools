classdef Mts
% class for storing a vector of timestamps (spike times for 1 cell)

    properties (GetAccess = public, SetAccess = private)
        t(:,1) double
        nspike(1,1) double
        tlim(1,2) double
    end

    methods (Access = public)
        function o = Mts(t)
        % constructor.  
        % o = Mts(t)
        % 
        % t - vector of timestamps in seconds
        %
        % o - Mts object
            o.t = t;
            o.nspike = numel(t);
            o.tlim = [min(t) max(t)];
        end

        function r = restrict(o,tr)
        % restricts s-matrix to times in 'tr'.
        % r = restrict(o,tr)
        %
        % tr - n x 2 start/end times (seconds)
        %
        % r - restricted s-matrix
            arguments (Input)
                o (1,1) Mts
                tr (:,2) double
            end

            n = size(tr,1);
            ns = numel(o.t);
            ix = false(ns,1);
            for i = 1:n
                ix = ix | (o.t>=tr(i,1) & o.t<=tr(i,2));
            end
            r = Mts(o.t(ix));
        end
    end

end
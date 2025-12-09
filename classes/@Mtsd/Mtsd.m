classdef Mtsd
    % class for timestamped data.  contains a time vector and voltage
    % values.  Time vector may be discontinuous, and voltage may be single
    % channel or multi-channel

    properties (GetAccess = public, SetAccess = private)
        t(:,1) {mustBeFloat} % time vector
        v {mustBeFloat} % voltage, nsample x nchannel
        fs double % sampling frequency
        tlim (1,2) double % start/end times
    end


    methods (Access = public)
        function obj = Mtsd(t,v)
        % constructor.  
        % obj = Mtsd(t,v)
        %
        % t - time vector in seconds
        % v - signal data (nsample x nchannel) 
        % t,v are forced to column arrangement, with the assumption that
        % there are more samples than channels.
        %
        % obj - Mtsd object
            if length(t) ~= length(v)
                error('length of time and signal must match')
            end
            obj.t = t;
            [nr,nc] = size(v);
            if nr>nc
                obj.v = v;
            else
                obj.v = v';
            end
            obj.fs = calcfs(obj);
            obj.tlim = [obj.t(1) obj.t(end)];
        end

        function o = plus(o1,o2)
        % concatenate 2 Mtsd's.
        % o = o1 + o2
        %
        % o1,o2 - Mtsd objects
        %
        % o - concatenated Mtsd
            o = Mtsd([o1.t; o2.t], [o1.v; o2.v]);
        end

        function c = getchan(o,ix)
        % Extract channel(s) from multichannel Mtsd
        % c = getchan(o,ix)
        %
        % ix - index(s) of channel(s) to extract. Can be numeric or
        %      logical.
        %
        % c - Mtsd with extracted channels
         
            if size(o.v,2) == 1
                c = o;
                warning('input has only 1 channel')
            else
                c = Mtsd(o.t,o.v(:,ix));
            end
        end

        function r = restrict(o,tr)
        % restricts 'o' to times in 'tr'
        % r = restrict(o,tr)
        % 
        % tr - n x 2 start/end times in seconds
        %
        % r - Mtsd of segment(s) in 'tr'
            arguments (Input)
                o (1,1) Mtsd
                tr (:,2) double
            end

            n = size(tr,1);
            ix = cell(1,n);
            for i = 1:n
                a = binsearch(o.t,tr(i,1));
                b = binsearch(o.t,tr(i,2));
                ix{i} = a:b;
            end
            ix = [ix{:}]';
            r = Mtsd(o.t(ix), o.v(ix,:));
        end

        function f = filter(o,flo,fhi,ord)
        % zero-phase butterworth filter.
        % f = filter(o,flo,fhi,ord)
        % 
        % flo,fhi - cutoff frequencies in Hz [0,inf]
        % ord -  optional filter order (default 8)
        %
        % f - filtered Mtsd
            arguments (Input)
                o (1,1) Mtsd
                flo (1,1) double {mustBeNonnegative}
                fhi (1,1) double {mustBeNonnegative}
                ord (1,1) double = 4
            end

            fny = o.fs * .5;
            hi = fhi/fny;
            lo = flo/fny;
            if lo == 0
                [b,a] = butter(ord,hi,'low');
            elseif hi == inf
                [b,a] = butter(ord,lo,'high');
            else
                [b,a] = butter(ord,[lo hi],'bandpass');
            end
            y = filtfilt(b,a,o.v);
            f = Mtsd(o.t,y);
        end

        function n = notch(o,f)
            % notch filter.
            % n = notch(f)
            % 
            % f - optional frequency to notch, default is 60 Hz
            %
            % n - filtered Mtsd
            arguments (Input)
                o (1,1) Mtsd
                f double = 60
            end
            if isscalar(f)
                fc = [f-1 f+1];
            else
                fc = f;
            end
            fn = o.fs/2;
            w = fc/fn;
            [b,a] = butter(4,w,'stop');
            y = filtfilt(b,a,o.v);
            n = Mtsd(o.t,y);
        end

        function d = downsample(o,fsnew)
        % downsamples data by decimation after applying anti-alias filter
        % d = downsample(o,fsnew)
        %
        % fsnew - desired sampling frequency in Hz
        % 
        % d - downsampled Mtsd.  New sampling frequency is limited by
        %     decimation and will be <= fsnew
            if fsnew >= o.fs
                d = o;
                warning('requested sampling rate is >= current sampling rate')
                return
            end
            fny = fsnew * .5;
            d = o.filter(0,fny);
            n = floor(o.fs/fsnew);
            d = Mtsd(d.t(1:n:end), d.v(1:n:end,:));
        end

        function montage(o,tp,scale,ax)
        % Plots trace(s) with optional timestamp patched.  
        % Multiple channels are spaced vertically.
        % montage(o,tp)
        %
        % tp - optional timestamps to patch (nx2 start/end times in seconds)
            n = size(o.v,2);
            if n > 1
                o.v = zscore(o.v);
                y = 0:20:20*(n-1);
                if nargin>=3
                    y = y.*scale;
                end
                o.v = bsxfun(@plus,o.v,y);
            end
            if nargin==4
                axes(ax)
                hold on
            else
                figure;
                tiledlayout(1,1,'padding','compact')
                nexttile;
            end
            plot(o.t,o.v,'k')
            xlabel('Time (sec)')
            if n > 1
                set(gca,'ytick',y,'yticklabel',1:n)
            end
            axis tight
            if nargin>1 && ~isempty(tp)
                tp = restrict_ts(tp,o.tlim);
                x = [tp fliplr(tp)]';
                y = ylim;
                y = [y;y];
                y = repmat(y(:),1,size(x,2));
                patch(x,y,'r','edgecolor','none','facealpha',.25);
            end
        end

        function [p,f,t] = spect(o,windur,flim,blur)
        % makes spectrogram image.
        % [p,f,t] = spect(o, windur, flim, blur)
        %
        % windur - duration of window in sec, default is 2.  
        %          Window is a Hann window with 50% overlap between 
        %          adjacent samples.  
        % flim - Frequency limits (Hz) for plot, default is [0 30]
        % blur - blur radius in standard deviations.  Default is [.5 1]
        %        Specify empty for no blur.
        %
        % p - freqency x time power matrix
        % f - frequency values
        % t - time bin values
            if size(o.v,2) > 1
                error('must be single channel')
            end
            if nargin<2 || isempty(windur)
                windur = 2;
            end
            if nargin<3 || isempty(flim)
                flim = [0 30];
            end
            if nargin<4
                blur = [.5 1];
            end
            
            win = round(windur*o.fs);
            w = hann(win);
            npt = max(4096,win);
            olap = round(win/2);
            [p,f,t] = spectrogram(o.v, w, olap, npt, o.fs);
            ix = f>=flim(1) & f<=flim(2);
            p = abs(p(ix,:)).^2 ./ npt^2;
            f = f(ix);
            if ~isempty(blur)
                p = imgaussfilt(p,blur);
            end
            t = t + o.t(1);
            tu = 'sec';
            if range(t) > 120
                t = t/60;
                tu = 'min';
            end
            figure;
            tiledlayout(1,1,'padding','compact')
            nexttile;
            imagesc(t,f,p)
            set(gca,'ydir','normal')
            xlabel(['Time (' tu ')'])
            ylabel('Frequency (Hz)')
        end

        function [av,t,se] = trigavg(o,tr,dur,noavg)
            arguments (Input)
                o (1,1) Mtsd
                tr (:,1) double
                dur (1,1) double
                noavg (1,1) logical = false
            end
            ix = tr > o.tlim(1)+dur/1000 & tr < o.tlim(2)-dur/1000;
            tr = tr(ix);
            n = numel(tr);
            dt = 1000/o.fs;
            nbin = round(dur/dt);
            a = zeros(n,2*nbin+1);
            for i = 1:n
                ix = binsearch(o.t,tr(i));
                a(i,:) = o.v(ix-nbin:ix+nbin)';
            end
            if noavg
                av = a;
            else
                av = mean(a);
            end
            if nargout>1
                t = (-nbin:nbin)*dt;
                if nargout>2
                    se = std(a) ./ sqrt(n);
                end
            end
        end
    end


    methods (Access = private)
        function srate = calcfs(o)
        % calculate sampling rate, warn of discontinuities
        % greater than 10% of dt (1/fs)
            dt = diff(o.t);
            md = median(dt);
            srate = 1/md;
            % round sampling rate if it is off by tiny amount
            d = abs(srate-round(srate));
            if d < 1e-4
                srate = round(srate);
            end
            thresh = md + md*.1;
            ix = dt > thresh;
            if any(ix)
                mx = max(dt(ix));
                u = 'sec';
                if mx < 1
                    mx = mx * 1000;
                    u = 'ms';
                elseif mx > 60
                    mx = mx / 60;
                    u = 'min';
                end
                s = 'Possible discontinuity in timestamps';
                s = sprintf('%s\nFound max gap of: %.1f %s\n',s,mx,u);
                warning(s);
            end
        end
    end
end


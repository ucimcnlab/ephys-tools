function hdr = intan_settings(fname)
%
% Read recording settings from intan settings xml file.  File name argument
% is only necessary if name is not 'settings.xml'.  Requires the file
% 'parseXML.m' to work.
% 
% Returns settings in the struct 'hdr':
%
%                                 fs: 30000
%                   AnalogScaleVolts: 1
%                         DSPEnabled: 1
%          DesiredDSPCutoffFreqHertz: 1
%     DesiredLower3dBCutoffFreqHertz: 1
%         DesiredLowerBandwidthHertz: 0.1000
%         DesiredUpperBandwidthHertz: 7500
%                              nchan: [1×4 struct]
%

validpfx = 'ABCD';
atts = {'AnalogScaleVolts','DSPEnabled','DesiredDSPCutoffFreqHertz',...
    'DesiredLower3dBCutoffFreqHertz','DesiredLowerBandwidthHertz','DesiredUpperBandwidthHertz'};

if nargin<1
    fname = 'settings.xml';
end
a = parseXML(fname);
c = {a.Attributes.Name};
ix = strcmp(c,'SampleRateHertz');
hdr.fs = str2double(a.Attributes(ix).Value);
chld = {a.Children.Name};
ix = strcmp(chld,'GeneralConfig');
b = a.Children(ix);
c = {b.Attributes.Name};
for i = 1:length(atts)
    ix = strcmp(c,atts{i});
    v = b.Attributes(ix).Value;
    hdr.(atts{i}) = str2num(lower(v));
end

nch = struct;
ig = 0;
si = find(strcmp(chld,'SignalGroup'));
fprintf('\nValid probe channels by shank:\n');
ntot = 0;
for i = 1:length(si)
    b = a.Children(si(i));
    pfx = {b.Attributes.Name};
    ix = strcmp(pfx,'Prefix');
    pfx = b.Attributes(ix).Value;
    ix = ismember(validpfx,pfx);
    if isscalar(pfx) && any(ix)
        ig = ig + 1;
        nch(ig).prefix = validpfx(ix);
        b = b.Children;
        ix = strcmp({b.Name},'Channel');
        b = b(ix);
        n = length(b);
        ok = false(1,n);
        num = nan(1,n);
        for j = 1:n
            att = {b(j).Attributes.Name};
            ix = strcmp(att,'NativeChannelName');
            name = b(j).Attributes(ix).Value;
            v = str2double(name(3:end));
            ix = strcmp(att,'Enabled');
            en = strcmp(b(j).Attributes(ix).Value,'True');
            if ~isnan(v) && en
                ok(j) = true;
                num(j) = v;
            end
        end
        nch(ig).prefix = pfx;
        nch(ig).chans = num;
        nch(ig).okchan = ok;
        fprintf('%s: %d\n',pfx,nnz(ok));
        ntot = ntot + nnz(ok);
    end
    1;
end
fprintf('Total: %d\n\n',ntot);
hdr.nchan = nch;
1;
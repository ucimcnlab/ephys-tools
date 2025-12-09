function data = Preprocessing(mouseName,garrFile,varargin)
% data = Preprocessing(mouseName,garrFile,intanFolder)
%
% Parses beahvioral data (from garrFile and intan) and corrects any time 
% difference between garrFile and intan.
% 
%   required inputs:
%   mouseName: name of mouse (string) for prefix of intan files
%   garrFile: name of file where garr mat data is stored (char)
%
%   optional name, value inputs:
%   'VRradius' - radius of ball (cm), default: 50
%   'samplerate' - sample rate of behavior data (Hz), default: 30
%   'nbin' - number of bins to divide track into, default: 100
%   'intanFolder' - folder where intan files are stored, 
%                   default: current folder
%   'doplot' - true/false for generating plots, default: true
%
%   Output:
%   data: structure containing behavior data.  Fields are:
%           dataLap: table containing all raw data from garr file
%    nonzeroVel_idx: indices of movement in raw data
%        dataBinned: structure with data binned by position
%                      posbins: 1 x nbin+1 bin edges used to segment track
%                     timebins: nlap x nbin matrix of time bins
%                               corresponding to position
%                    occupancy: nlap x nbin matrix of time spent in each
%                               bin (s)
%                     velocity: nlap x nbin matrix of velocity in each bin
%                               (cm/s)
%              date: date extracted from name of garr file
%       intanoffset: offset in samples of intan relative to garr file
%       time_offset: time offset in seconds
%            reward: nx2 matrix with time,position of rewards
%          garrfile: file name data came from.
%
% RL 11/2024 - original
% ME 3/2025 - fixed issue with non-sequential bin edges, 
%             changed representation of binned data,



% Parameters
par = inputParser;
addParameter(par, 'VRradius',50,@isnumeric);
addParameter(par, 'samplerate',30,@isnumeric);
addParameter(par, 'nbin',100,@isnumeric);
addParameter(par, 'intanFolder',pwd, @isfolder);
addParameter(par, 'doplot',true, @(x)isscalar(x)&&islogical(x));
parse(par, varargin{:});
par = par.Results;

VRradius = par.VRradius;
samplerate = par.samplerate;
nBin = par.nbin;
intanFolder = par.intanFolder;
doplot = par.doplot;
binsize = 1/nBin;

% check if data files exist
if ~isstring(mouseName)
    mouseName = string(mouseName);
end
if ~isfile(garrFile)
    error('garr file not found')
end
fn_digin = fullfile(intanFolder,mouseName+'-digIn.npy');
if ~isfile(fn_digin)
    error('file not found: %s',fn_digin)
end
fn_digints = fullfile(intanFolder,mouseName+'-digInts.npy');
if ~isfile(fn_digints)
    error('file not found: %s',fn_digints)
end

% load data
load(garrFile, 'garr', 'adata');
if garr(end,1) == 0
    garr(end,:) = [];
end
s = strsplit(garrFile,'_');
ix = find(contains(s,digitsPattern),1,'first');
fdate = s{ix};
digIn = readNPY(fn_digin);
digInts = readNPY(fn_digints);
intants = loadintants(fn_digin);

mkdir Preprocess_figures % To store saved figures

% TODO 
% not sure about this, seems to find radius for each env, but radius
% should be fixed since it's the ball?  Ignore for now
environment = garr(:,7);
[envs,~,ic] = unique(environment);
% for ii = 1:length(envs)
%     [~,VRradius(ii)] = findEnvironment(envs(ii));
% end
%vrrad_eachbin = VRradius(ic);  % radius of environment at each time point, seems to be very accurate even at env transitions

nsample = size(garr,1);
pos_raw = garr(:,5); % position in laps (based on VR)
poscm_raw = pos_raw*VRradius*2*pi;
pos = mod(pos_raw,1);
poscm = pos*VRradius*2*pi;
posBin = floor(pos*nBin) + 1;
samples = [0:nsample-1]';
vel = [0; diff(poscm_raw)] * samplerate; % Velocity (cm/s)
vel_smoothed = smoothvel(vel);
ts = samples / samplerate; % Timestamp (s)

reward = [0; diff(garr(:,3)) < 0]; % Reward
licks = garr(:,4); % Lick
objectChanges = garr(:,6);


% Redundant variables (for now)
% pos_enc = garr(:,1); % position in revolutions of wheel (based on rotary encoder)
% pos_incm_enc = pos_enc.*12.*pi; % position in cm (based on rotary encoder) -- in movie mode this will be different from pos
% inst_vel = garr(:,2);  %cm/s (based on VR)
% vel_enc = [0;diff(pos_incm_enc)];   %cm/sample (based on encoder)
% vel_enc = vel_enc.*samplerate; %cm/s (based on encoder)

% Split dataset into laps
lapStart_idx = [1; find(diff(pos) < -0.95) + 1]; % Threshold difference at 95% of track distance
lapEnd_idx = find(diff(pos) < -0.95);
nLap = length(lapEnd_idx);
lapStart_idx = lapStart_idx(1:nLap);

lap = ones(nsample,1);
for lap_no = 1:nLap
    lap(lapStart_idx(lap_no):lapEnd_idx(lap_no)) = lap_no;
end
lap(lapEnd_idx(end)+1:end) = NaN;

dataLap_all = [lap posBin ts pos poscm vel vel_smoothed reward licks objectChanges environment];

% Need movement data only (zero velocity)
posChange = (diff(dataLap_all(:,4)) ~= 0);
nonzeroVel_idx = logical([posChange; 0]);
dataLap = dataLap_all(nonzeroVel_idx,:);

% Align timestamps with Intan
dt = diff(dataLap(:,3));
dint = diff(intants);
[c,lag] = xcorr(dint,dt);
[mx,ix] = max(c);
offset = lag(ix);
if offset > 0 % Intan starts first
    time_offset = intants(abs(offset));
    offset_shift = find(intants > time_offset,1,'first');
    intants = intants(offset_shift:end);
else          % Smoothwalk starts first
    time_offset = dataLap(abs(offset),3);
    offset_shift = find(dataLap(:,3) > time_offset,1,'first');
    dataLap = dataLap(offset_shift:end,:);
end

% Plot for timestamp offset between behavior and intan
if doplot
    f = figure;
    tp = tiledlayout(2,1);
    nexttile;
    plot(dt);
    hold on
    plot(dint);
    toff = offset/samplerate;
    if abs(toff) > 120
        toff = toff/60;
        lab = 'min';
    else
        lab = 'sec';
    end
    title(sprintf('Intan timestamp offset: %.0f (%.1f %s)',offset,toff,lab))
    xlabel("Timestamps"); ylabel("Magnitude of Difference (s)");
    legend('VR','intan')
    nexttile;
    plot(lag,c);
    hold on
    h = plot([0 0],[0 mx*1.05],'--','color',[.7 .7 .7],'linewidth',1.25);
    h.ZData = [-1 -1];
    axis tight
    title('xcorr of intan and VR')
    xlabel('intan difference from VR (samples)');
    tp.TileSpacing = 'compact'; tp.Padding = 'compact';
    f.Position = f.Position + [0 -100 -100 100];
    saveas(gcf,"Preprocess_figures/" + mouseName +  "_Timestamp_Alignment.png")
end

dataLap(:,3) = intants;

% Remove fake laps (due to jittering of mice position around lap transition point)
ibad = diff(dataLap(:,2)) > 90;
fakeLaps = unique(dataLap(ibad,1));
fakeLaps = [dataLap(1,1); fakeLaps]; % Also remove the first incomplete lap
fakeLaps_idx = zeros(size(dataLap,1),1);
fakeLapNo_correction = dataLap(:,1);
for i = 1:length(fakeLaps)
    temp_idx = find(dataLap(:,1) == fakeLaps(i));
    fakeLaps_idx(temp_idx) = 1;
    fakeLapNo_correction(temp_idx(1):end) = fakeLapNo_correction(temp_idx(1):end) - 1; % Fix lap no.
end
if dataLap(1,1) > 1 % Start lap is not labeled as 1
    fakeLapNo_correction = fakeLapNo_correction - dataLap(1,1) + 1;
end
fakeLaps_idx(isnan(dataLap(:,1))) = 1; % Remove final incomplete lap
nLap = nLap - (dataLap(1,1) - 1) - length(fakeLaps);
dataLap(:,1) = fakeLapNo_correction;
dataLap(logical(fakeLaps_idx),:) = []; % Remove fake laps + first and final incomplete lap from behaviour data

% Bin velocity by position bin per lap
% mike new code to avoid time bin errors
dataBinned.posbins = 0:binsize:1;
dataBinned.timebins = zeros(nLap,nBin+1);
dataBinned.occupancy = zeros(nLap,nBin);
dataBinned.velocity = zeros(nLap,nBin);
for i = 1:nLap
    ix = dataLap(:,1) == i;
    posind = dataLap(ix,2);
    t_ = dataLap(ix,3);
    v_ = dataLap(ix,7);
    ed = accumarray(posind,t_,[],@min);
    v = accumarray(posind,v_,[],@(x)mean(x,'omitnan'));
    ii = find(ed>0);
    if length(ii)<nBin
        ed = interp1(ii,ed(ii),1:nBin)';
        v = interp1(ii,v(ii),1:nBin);
    end
    ed = [ed; max(t_)];
    dataBinned.timebins(i,:) = ed;
    dataBinned.occupancy(i,:) = diff(ed);
    dataBinned.velocity(i,:) = v;
end

% royston original code for binning
% for i = 1:size(dataBinned,1)
%     lap_no = floor((i-1)/nBin) + 1;
%     dataBinned(i,1) = lap_no; % Lap No.
%     bin_no = mod(i,nBin);
%     if mod(i,nBin) == 0
%         bin_no = nBin;
%     end
%     dataBinned(i,2) = bin_no; % Bin No.
%     temp_array = dataLap((dataLap(:,1) == lap_no) & (dataLap(:,2) == bin_no),:);
%     if ~isempty(temp_array)
%         dataBinned(i,3) = temp_array(1,3); % Start time
%         dataBinned(i,4) = temp_array(end,3); % End time
%         dataBinned(i,5) = temp_array(end,3) - temp_array(1,3) + (1/samplerate); % Duration
%         dataBinned(i,6) = mean(temp_array(:,7),'omitnan'); % Use smoothed velocity
%     end
% end


% Extract reward signals from intan
digInts_arr = digInts(find(digInts >= dataLap(1,3),1,'first'):find(digInts <= dataLap(end,3),1,'last'));
digIn_arr = digIn(2,find(digInts >= dataLap(1,3),1,'first'):find(digInts <= dataLap(end,3),1,'last'));
digIn_arr_ = [0 (abs(diff(int8(digIn_arr))))];
intan_pos = interp1(dataLap(:,3),dataLap(:,5),digInts_arr);

reward_ts = dataLap(dataLap(:,8) == 1,3);
reward_pos = dataLap(dataLap(:,8) == 1,5) / (VRradius*2*pi);
intan_reward_ts = digInts_arr(digIn_arr_ == 1);
intan_reward_pos = intan_pos(digIn_arr_ == 1) / (VRradius*2*pi);

% plot reward position
if doplot
    figure;
    subplot(2,1,1)
    scatter(reward_ts,reward_pos,'bx'); yline(0.83,'k--');
    title("Smoothwalk"); xlabel("Time (s)"); ylabel("Reward Position");
    ylim([0 1]);
    subplot(2,1,2)
    scatter(intan_reward_ts,intan_reward_pos,'bx'); yline(0.83,'k--');
    title("Intan"); xlabel("Time (s)"); ylabel("Reward Position");
    ylim([0 1]);
    saveas(gcf,"Preprocess_figures/" + mouseName +  "_Lap_Reward_Position.png")
end


% put all relevant data in structure and save
vars = {'lapnum' 'posBin' 'time' 'pos' 'poscm' 'vel' 'vel_smoothed'...
    'reward' 'licks' 'objectChanges' 'environment'};
data.dataLap = array2table(dataLap_all,'variablenames',vars);
data.nonzeroVel_idx = nonzeroVel_idx;
data.dataBinned = dataBinned;
data.date = fdate;
data.intanoffset = offset;
data.time_offset = time_offset;
data.reward = [intan_reward_ts intan_reward_pos];
data.garrfile = garrFile;
fn = mouseName + '_behavior_' + fdate + '.mat';
save(fn,"data",'-v7.3');

if doplot
    % Plot behavioural data
    f = figure;
    tp = tiledlayout(2,2);
    % First 10 laps
    start_lap = 1;
    end_lap = 10;
    plot_start_idx = find(dataLap(:,1) == start_lap,1,'first');
    plot_end_idx = find(dataLap(:,1) == end_lap,1,'last');
    nexttile;
    yyaxis left
    plot(dataLap(plot_start_idx:plot_end_idx,3),dataLap(plot_start_idx:plot_end_idx,5))
    title("Behavior - first 10 laps");
    xlabel('Time (s)'); ylabel('Distance (cm)');
    hold on
    yline(0.83*VRradius*2*pi,'--');
    yyaxis right
    digInts_arr = digInts(find(digInts >= dataLap(plot_start_idx,3),1,'first'):find(digInts <= dataLap(plot_end_idx,3),1,'last'));
    digIn_arr = digIn(2,find(digInts >= dataLap(plot_start_idx,3),1,'first'):find(digInts <= dataLap(plot_end_idx,3),1,'last'));
    digIn_arr_ = [0 (abs(diff(int8(digIn_arr))))];
    plot(digInts_arr,digIn_arr_,'r')
    ylabel('Reward');
    nexttile(3);
    plot(dataLap(plot_start_idx:plot_end_idx,3),dataLap(plot_start_idx:plot_end_idx,6))
    title("Velocity (cm/s)");
    xlabel('Time (s)'); ylabel('Velocity (cm/s)');
    ylim([0,120])
    hold on
    % Last 10 laps
    start_lap = dataLap(end,1) - 9;
    end_lap = dataLap(end,1);
    plot_start_idx = find(dataLap(:,1) == start_lap,1,'first');
    plot_end_idx = find(dataLap(:,1) == end_lap,1,'last');
    nexttile(2);
    yyaxis left
    plot(dataLap(plot_start_idx:plot_end_idx,3),dataLap(plot_start_idx:plot_end_idx,5))
    title("Behavior - last 10 laps");
    xlabel('Time (s)'); ylabel('Distance (cm)');
    hold on
    yline(0.83*VRradius*2*pi,'--');
    yyaxis right
    digInts_arr = digInts(find(digInts >= dataLap(plot_start_idx,3),1,'first'):find(digInts <= dataLap(plot_end_idx,3),1,'last'));
    digIn_arr = digIn(2,find(digInts >= dataLap(plot_start_idx,3),1,'first'):find(digInts <= dataLap(plot_end_idx,3),1,'last'));
    digIn_arr_ = [0 (abs(diff(int8(digIn_arr))))];
    plot(digInts_arr,digIn_arr_,'r')
    ylabel('Reward');
    nexttile(4)
    plot(dataLap(plot_start_idx:plot_end_idx,3),dataLap(plot_start_idx:plot_end_idx,6))
    title("Velocity (cm/s)");
    xlabel('Time (s)'); ylabel('Velocity (cm/s)');
    ylim([0,120])
    tp.TileSpacing = 'compact'; tp.Padding = 'compact';
    f.Position = f.Position + [0 -100 400 100];
    saveas(gcf,"Preprocess_figures/" + mouseName +  "_Distance_Velocity.png")

    % Plot position-binned data
    f = figure;
    imagesc(dataBinned.velocity); cbar = colorbar; cbar.Label.String = 'Velocity (cm/s)';
    f.Position = f.Position + [0 -100 -100 100];
    title("Velocity (cm/s)");
    xlabel('Position (Bins)'); ylabel('Lap');
    saveas(gcf,"Preprocess_figures/" + mouseName +  "_Velocity_Binned.png")

    figure;
    xdat = (1:nBin)*binsize*VRradius*2*pi;
    ydat = mean(dataBinned.velocity,'omitnan');
    y_std = std(dataBinned.velocity,'omitnan');
    plot(xdat,ydat,'k-','LineWidth', 2);
    hold on
    curve1 = ydat + y_std;
    curve2 = ydat - y_std;
    x2 = [xdat, fliplr(xdat)];
    inBetween = [curve1, fliplr(curve2)];
    fill(x2, inBetween, 'k','FaceAlpha',0.3);
    title(mouseName);
    xlabel('Position (cm)'); ylabel('Velocity (cm/s)');
    xlim tight
    saveas(gcf,"Preprocess_figures/" + mouseName +  "_Velocity.png")
end

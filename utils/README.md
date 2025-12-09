# ephys-utils

### File importing  
**eeg = npy2mtsd(fname, fs, chan, chunk)**  - Reads a numpy LFP file and makes a Mtsd object.  This function requires the 'npy-matlab' package as well as 'readnpychunk' and 'readnpychan'.   
```
fname - numpy file name containing LFP data

Optional name/value pairs
'fs' - sampling frequency (Hz). If not given, tries to find timestamp file.
'channel' - optional single channel to extract, default is all
'chunk' - optional start/stop indices for import, default is whole file
```
<br>

**a = readnpychunk(fname, chunk)**  - Reads a subset of a npy LFP file and returns the raw data as 2D matrix.  
```
Read a portion of a numpy file.  Reads all channels.
Works for 2D single or double matrices

fname - file name
ix - 2 element vector specifying start/end indices (samples), defaults
     to whole recording if not given / empty.
```
<br>  

**a = readnpychan(fname, chan)**  - Reads a single channel of npy LFP file and returns a vector.  This is slow.  
```
Read a single channel from a numpy file
Works for 2D single or double matrices

fname - Name of numpy file
ch - channel number to extract
```
<br>  

**[eeg,par] = intan2mtsd(fname, varargin)**  - Reads raw data from the intan .bin file, converts to mV, and makes a Mtsd object 'eeg'.  Also reads the 'settings.xml' file to determine sampling rate various other recording parameters and returns them in the 'par' structure.
```
required input:
fnbin - file name of binary file

optional name,value args:
'channels' - which channels to plot, default: [] (empty=all)
'epoch' - time range to plot in sec: [start end], default: [] (all)
'doread' - true/false if you really want to read the file.  Otherwise,
           information about the file is returned. default: true
'settingsfile' - name of xml settings file, default: 'settings.xml'

Output
eeg - data as a Mtsd object
par - structure with file information and parameters. Fields are:
         binaryfile - name of binary file
      setttingsfile - name of settings file
      validchannels - all valid channels on the probe
           channels - requested channels to read
              epoch - time interval specified for reading
                 fs - sampling frequency of the data
              scale - multiplier used to convert data to volts
           nsamprec - number of data samples in the file
             recdur - duration of recording in sec
             doread - true/false indicating if read was requested
```
<br>  

**eeg = csc2mtsd(fname, ts)**  - imports a neuralynx CSC file as a Mtsd object.  
```
 Read neuralynx csc file (LFP data) and return Mtsd object.  

 fname - name of csc file
 ts - optional start/end timestamps (seconds) to extract. 
      default: [] (whole file)
```
<br>  

**vt = vt2mtsd(varargin)**  - import neuralynx tracker coordinates as a Mtsd object.
```
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


vt = vt2mtsd
vt = 

  Mtsd with properties:

       t: [11132×1 double]
       v: [11132×2 double]
      fs: 30.0210
    tlim: [1.0669e+06 1.0673e+06]

```

### Load behavior data from ball rig  

**dat = Preprocessing(mouseName, garrFile, varargin)**  - Parses behavioral data from behavior and ephys files and returns the structure 'dat'.  Also checks and corrects any timestamp offset between the two data files.  
```
   required inputs:
   mouseName: name of mouse (string) for prefix of intan files
   garrFile: name of file where garr mat data is stored (char)

   optional name, value inputs:
   'VRradius' - radius of ball (cm), default: 50
   'samplerate' - sample rate of behavior data (Hz), default: 30
   'nbin' - number of bins to divide track into, default: 100
   'intanFolder' - folder where intan files are stored, 
                   default: current folder
   'doplot' - true/false for generating plots, default: true

   Output:
   dat: structure containing behavior data.  Fields are:
           dataLap: table containing all raw data from garr file
    nonzeroVel_idx: indices of movement in raw data
        dataBinned: structure with data binned by position
                      posbins: 1 x nbin+1 bin edges used to segment track
                     timebins: nlap x nbin matrix of time bins
                               corresponding to position
                    occupancy: nlap x nbin matrix of time spent in each
                               bin (s)
                     velocity: nlap x nbin matrix of velocity in each bin
                               (cm/s)
              date: date extracted from name of garr file
       intanoffset: offset in samples of intan relative to garr file
       time_offset: time offset in seconds
            reward: nx2 matrix with time,position of rewards
          garrfile: file name data came from.
```
The script also generates several plots by default: the average and individual lap running speed, the timestamp offset between the behavior and ephys datasets, the reward position across laps, and a comparison of the position, velocity and reward loctation for the first 10 and last 10 laps.  
<img src="/images/behav_vel_trials.png" width=30%><img src="/images/behav_vel_avg.png" width=40%>
<img src="/images/behav_timstamp_offset.png" width=30%>
<img src="/images/behav_first_last_laps.png" width=50%>
<br><br>

### Load Kilosort spike data  

**dat = load_spike_data(varargin)**  - Loads spike data from Kilosort and optionally classifies clusters as pyramidal or interneuron based on firing rate and waveform.  
```
optional name,value arguments:
'fn_info' - name of kilosort cluster info file. default: 'cluster_info.tsv'
'fn_clu' - name of kilosort spike cluster file. default: 'spike_clusters.npy'
'fn_times' - name of kilosort spike times file. default: 'spike_times.npy'
'fs' - sampling frequency in Hz. default: 30,000
'groups' - which group labels to extract.  Can be a single group or a 
           cell array. default: 'good'
'epoch' - limit data extraction to time interval as a [1 2] vector of
          start/end times in seconds.  Default is to use global min/max
          of all spiketimes pooled across cells.
'doclassify' - true/false if you want to classify clusters as pyramidal
               or interneuron based on peak-trough duration and firing
               rate.  Requires bombcell waveform duration file.
'class_fr' - firing rate threshold for classification in Hz. default: 8
'class_dur' - waveform peak to trough duration for classification in
              microseconds. default: 400
'fn_wav' - name of bombcell waveform duration file.
           default: 'cluster_waveform_duration.tsv'

Outputs structure 'dat' with fields:
        clu: cluster number
      depth: position of cell along electrode (microns)
         ch: channel of max amplitude
 spiketimes: spiketimes in seconds
         fr: average firing rate (Hz)
   celltype: 'pyramidal' or 'interneuron' (if requested)
```
<br>

**hdr = intan_settings(fname)**  - Read recording settings from intan configuration file.  
```
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
```
<br>

**S = parseXML(fname)**  - General purpose XML reader that formats data in a struct.  
```
S = parseXML('settings.xml')

  struct with fields:

          Name: 'IntanRHX'
    Attributes: [1×3 struct]
          Data: ''
      Children: [1×17 struct]
```
<br>

### Analysis  

**maps = ratemap1d(tbin,vel,spikes,varargin)**   - calculates 1D firing rate maps for a set of cells.  Designed for behavior on a linear track or treadmill.  
```
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
```
<br>

**ax = ratemap_figs(rmap, varargin)**  - plots 1D ratemaps over trials as images with dimensions ntrial x nbin and color for firing rate in Hz.
```
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
```
Example ratemap  
<img src="/images/ratemap.png" width=70%>
<br>


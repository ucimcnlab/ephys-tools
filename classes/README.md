# ephys-classes
Classes for working with spike and LFP data

Classes for working with ephys data adapted from old McNaughton lab versions.

## Mts 
a vector of timestamps, e.g., spike times.
Create an Mts object:
```
ts = [1.1; 2.3; 3.2; 4.7];
mts_obj = Mts(ts);
```

#### Properties
**nspike** - the number of spikes  
**tlim** - min and max spike times  

#### Methods
**mts_obj = restrict(ts)** - extracts spikes between specified start/end times (sec) in 'ts' [1 x 2] and returns a new Mts.  
Restrict to times between 2 and 4 seconds:
```
mts_obj_res = mts_obj.restrict([2 4]);
```

## Msmat
A cell array of Mts objects, called an S-matrix, for working with the spike times from a number of cells.  
Construct a Msmat from raw spike times:  
```
spikes = cell(1,10);
spikes = cellfun(@(x) rand(1,100)*10, spikes, 'UniformOutput', false);
msmat_obj = Msmat(spikes)
msmat_obj = 

  Msmat with properties:

     cells: {1×10 cell}
     ncell: 10
      tlim: [0.0052 9.9949]
    nspike: [100 100 100 100 100 100 100 100 100 100]
```
Construct a Msmat from a cell array of Mts objects:  
```
mts_objects = {mts_obj1, mts_obj2, mts_obj3};
msmat_obj = Msmat(mts_objects);
```

#### Properties  
**cells** - the cell array of Mts objects  
**ncell** - the number of cells  
**tlim** - the min and max spike times of the Smatrix  
**nspike** - number of spikes for each cell returned as a vector

#### Methods  
**smat = restrict(ts)** - extracts spikes between start/end timestamps 'ts' (1 x 2) for each cell and returns a new Msmat object.  
```
smat_res = msmat_obj.restrict([2 4]);
```
**raster(tp)** - makes a raster plot of the Smatrix and optionally adds start/end timestamp patches given in 'tp' (nx2)  
```
msmat_obj = Msmat(spikes);
tp = [2 4; 6 7];
msmat_obj.raster(tp)
```
<img src="/images/msmat_raster.png" width=40%>


## Mqmat  
Matrix of binned spike counts, called a Q-matrix (nbin x ncell)  
Construct a Qmatrix from a Smatrix (or a cell array of raw spike times) and the desired bin size in seconds:  
```
>> qmat_obj = Mqmat(msmat_obj, 0.1)

qmat_obj = 

  Mqmat with properties:

    tstart: 0.0052
        dt: 0.1000
      data: [99×10 double]
      tlim: [0.0052 9.9052]
     ncell: 10
```
Without specifying start/end times, the default is to use the min/max spike times in the Smatrix.  Alternatively, you can specify start/end times:
```
qmat_obj = 

  Mqmat with properties:

    tstart: 1
        dt: 0.1000
      data: [30×10 double]
      tlim: [1 4]
     ncell: 10
```
You can also contstruct a Qmatrix from spike times that are already binned into a nbin x ncell matrix.  To do this, you must specify the start time and bin size in seconds:
```
>> spikecounts = round(rand(100,10)*10);
>> tstart = 0;
>> binsize = 0.1;
>> qmat_obj = Mqmat(tstart, binsize, spikecounts)

qmat_obj = 

  Mqmat with properties:

    tstart: 0
        dt: 0.1000
      data: [100×10 double]
      tlim: [0 10]
     ncell: 10
```

#### Properties  
**tstart** - start time in seconds  
**dt** - bin size in seconds  
**data** - binned spike counts (nbin x ncell).  **Note that data matrix is sparse, use full() to convert to double**  
**tlim** - start/end times for q-matrix  
**ncell** - number of cells  

#### Methods
**ed = binedges()** - returns the bin edge times in seconds  
**cen = bincenters()** - returns the times of bin centers in seconds  
**Qobj = restrict(ts)** - limits the Qmatrix to start/end times given in 'ts' (1 x 2) and returns a new Qmatrix  
**image(tp)** - makes an image of the Qmatrix with optional patches for start/end timestamps in 'tp' (npatch x 2):  
```
>> qmat_obj = Mqmat(msmat_obj, 0.1)

qmat_obj = 

  Mqmat with properties:

    tstart: 0.0052
        dt: 0.1000
      data: [99×10 double]
      tlim: [0.0052 9.9052]
     ncell: 10

>> tp = [2 4; 6 7];
>> qmat_obj.image(tp)
```
<img src="/images/qmatrix.png" width=60%>

**[av,t,se] = trigavg(tr,dur,noavg)** - make a triggered average of a Qmatrix  
Inputs  
tr - event times in seconds  
dur - duration to include pre/post trigger  
noavg - if true, do not calculate the average and return all events in a ncell x nbin x nevent matrix  
Outputs  
av - the triggered average as a ncell x nbin matrix  
t - time axis  
se - standard error of average  

Example triggered average with fake data:
```
% make some spike times
>> spikes = cell(1,10);
>> spikes = cellfun(@(x) rand(1,1000)*100, spikes, 'UniformOutput', false);
>> q = Mqmat(spikes,0.1)

q = 

  Mqmat with properties:

    tstart: 0.0070
        dt: 0.1000
      data: [999×10 double]
      tlim: [0.0070 99.9070]
     ncell: 10

% trigger times in seconds
>> tr = 5:5:95;

% increase spike counts in the second bin after trigger times
>> qdat = q.data;
>> bump = round(rand(length(tr),q.ncell)*3);
>> ix = 52:50:999;
>> qdat(ix,:) = qdat(ix,:) + bump;
>> qnew = Mqmat(q.tstart, q.dt, qdat);

% calculate and plot the triggered average
>> dur = 1;
>> [av,t,se] = qnew.trigavg(tr,dur);
>> figure; imagesc(t, [], av)
```
<img src="/images/trigavg.png" width=40%>

# Mtsd  
Class for LFP data that packages time and voltage vectors and provides convenient methods.  
Construct a Mtsd from time and voltage data.  Time and voltage must have the same length but voltage may contain multiple channels (nsample x nchannel):  
```
>> whos t v
  Name         Size             Bytes  Class     Attributes

  t         5001x1              40008  double              
  v         5001x12            240048  single              

>> eeg = Mtsd(t,v)

eeg = 

  Mtsd with properties:

       t: [5001×1 double]
       v: [5001×12 single]
      fs: 1000
    tlim: [2.6600e+03 2.6650e+03]
```

#### Properties  
**t** - time vector  
**v** - voltage data, [nsample x nchannel]  
**fs** - sampling frequency in Hz  
**tlim** - start/end times of data, [1 x 2]  

#### Methods  
**eeg3 = eeg1 + eeg2**  - concatenate 2 Mtsd objects with '+' operator  
**eeg = getchan(ch)**  - extracts subset of channels from multi-channel data and returns a new Mtsd, 'ch' may be boolean or vector of integers  
**eeg = restrict(ts)**  - extracts segment(s) of LFP based on start/end timestamps (seconds) in 'ts' [n x 2]  
**eeg = filter(flo, fhi, ord)** - zero-phase butterworth filter.  'flo' and 'fhi' are the low and high cutoff frequencies in Hz [1 x 2].  To make a low pass filter, specify 'flo' as 0, to make a high pass, specify 'fhi' as inf.  'ord' is optional order of filter (defualt 8).  
**eeg = notch(f)**  - notch filter for given frequency 'f' in Hz (default 60).  
**eeg = downsample(fsnew)**  - downsamples data by decimation after applying anti-alias filter.  'fsnew' is specified in Hz.  The new sampling rate is limited by decimation so it will be <= fsnew.  
**montage(tp, scale, ax)**  - plot the data and try to make multi-channel data look nice.  Note that data is zscored for convenient plotting.  'tp' is optional start/end timestamps to patch on top of the data [n x 2].  'scale' optionally controls the spacing between channels in multi-channel data, values <1 reduces spacing and values >1 increases spacing.  'ax' is optional axes to plot on.  
Example default montage of multiple LFP channels spanning hippocampus:  
```
>> eeg2.montage
```
<img src="/images/montage.png" width=60%>

**[p,f,t] = spect(windur, flim, blur)**  - uses the matlab function 'spectrogram' to image the time-frequency decomposition of a single channel Mtsd.  A sliding window is used to divide the signal into shorter segments before applying a FFT.  'windur' specifies the sliding window duration in seconds (default 2).  The window is a Hann window with 50% overlap between samples.  'flim' specifies the frequency limits in Hz of the resulting plot (default [0 30]).  'blur' specifies the standard deviations of the gaussian blur applied to the frequency and time dimensions [1 x 2].  Default is a mild blur of [.5 1].  For no blur, enter an empty matrix [].  Returns the power spectrum 'p' and the vectors for the frequency and time axes.  
Example spectrogram:  
```
>> eeg2 = eeg.getchan(6)

eeg2 = 

  Mtsd with properties:

       t: [5001×1 double]
       v: [5001×1 single]
      fs: 1000
    tlim: [4809 4814]

>> eeg2.spect(.15, [80 200], [.5 .7]);
```
<img src="/images/spect.png" width=40%>

**[av,t,se] = trigavg(tr, dur, noavg)**  - performs a triggered average of the LFP based on event times in the vector 'tr' (sec).  The length of the LFP segment to include around event times is given by 'dur' (sec).  'noavg' is a boolean that returns all LFP segments instead of the average if true.
Returns the average 'av' of the LFP segments [1 x nsample], or the matrix of segments [nsegment x nsample] if 'noavg' is true.  't' is the time axis [1 x nsample].  'se' is the standard error of the average [1 x nsample].


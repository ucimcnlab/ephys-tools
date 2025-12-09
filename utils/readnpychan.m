function a = readnpychan(fname,ch)
%
% Read a single channel from a numpy file
% Works for 2D single or double matrices
%
% fname - Name of numpy file
% ch - channel number to extract
%

[shape, dtype, fort, ~, hlen, ~] = readNPYheader(fname);

if strcmp(dtype,'double')
    nbyte = 8;
else
    nbyte = 4;
end

if shape(1) > shape(2) % data is in columns (assume more samples than channels)
    nch = shape(2);
    nread = shape(1);
    if fort
        rskip = 0;
        bskip = (ch-1)*nread*nbyte;
    else
        rskip = (nch-1)*nbyte;
        bskip = (ch-1)*nbyte;
    end
else % data is in rows
    nread = shape(2);
    nch = shape(1);
    if fort
        rskip = (nch-1)*nbyte;
        bskip = (ch-1)*nbyte;
    else
        rskip = 0;
        bskip = (ch-1)*nread*nbyte;
    end
end

fid = fopen(fname,'rb');
fseek(fid,double(hlen)+bskip,'bof');
a = fread(fid,nread,[dtype '=>' dtype],rskip);
fclose(fid);
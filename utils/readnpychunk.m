function a = readnpychunk(fname,ix)
%
% Read a portion of a numpy LFP file.  Reads all channels.
% Works for 2D single or double matrices
%
% fname - file name
% ix - 2 element vector specifying start/end indices (samples), defaults
%      to whole recording if not given / empty.
%

[shape, dtype, fort, ~, hlen, ~] = readNPYheader(fname);

if isscalar(shape)
    shape = [shape 1];
end

if strcmp(dtype,'double')
    nbyte = 8;
else
    nbyte = 4;
end

if nargin < 2 || isempty(ix)
    na = 1;
    nb = max(shape);
else
    na = ix(1);
    nb = ix(2);
end

nb = min(nb, max(shape));
n = nb-na+1;
if shape(1) > shape(2)
    nch = shape(2);
    nsamp = shape(1);
    rshape = [n nch];
    if fort
        nrem = nsamp-nb;
        bskip = (na-1)*nbyte;
        rskip = (nrem+na-1)*nbyte;
        fmt = sprintf('%d*%s=>%s',n,dtype,dtype);
    else
        bskip = nch*(na-1)*nbyte;
        rskip = 0;
        fmt = [dtype '=>' dtype];
    end
else
    nch = shape(1);
    nsamp = shape(2);
    rshape = [nch n];
    if fort
        bskip = nch*(na-1)*nbyte;
        rskip = 0;
        fmt = [dtype '=>' dtype];
    else
        nrem = nsamp-nb;
        bskip = (na-1)*nbyte;
        rskip = (nrem+na-1)*nbyte;
        fmt = sprintf('%d*%s=>%s',n,dtype,dtype);
    end
end
ntot = n*nch;

fid = fopen(fname,'rb');
fseek(fid,double(hlen)+bskip,'bof');
a = fread(fid,ntot,fmt,rskip);
fclose(fid);

if length(shape)>1 && ~fort
    a = reshape(a, rshape(end:-1:1));
    a = permute(a, [length(rshape):-1:1]);
elseif length(shape)>1
    a = reshape(a, rshape);
end
1;
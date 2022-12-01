function [signal,removed_channels] = rm_badchannels(signal,corr_threshold,noise_threshold,window_len,max_broken_time,num_samples,subset_size)

% Remove channels with abnormal data from a continuous data set.
% Signal = clean_channels(Signal,CorrelationThreshold,WindowLength,MaxBrokenTime,NumSamples,SubsetSize,UseGPU)
%
% This is an automated artifact rejection function which ensures that the data contains no channels
% that record only noise for extended periods of time. If channels with control signals are
% contained in the data these are usually also removed. The criterion is based on correlation: if a
% channel has lower correlation to its robust estimate (based on other channels) than a given threshold
% for a minimum period of time (or percentage of the recording), it will be removed.
%
% In:
%   Signal          : Continuous data set, assumed to be appropriately high-passed (e.g. >0.5Hz or
%                     with a 0.5Hz - 2.0Hz transition band).
%
%   CorrelationThreshold : Correlation threshold. If a channel is correlated at less than this value
%                          to its robust estimate (based on other channels), it is considered abnormal in
%                          the given time window. Default: 0.85.
%                     
%   LineNoiseThreshold : If a channel has more line noise relative to its signal than this value, in
%                        standard deviations from the channel population mean, it is considered abnormal.
%                        Default: 4.
%
%
%   The following are detail parameters that usually do not have to be tuned. If you cannot get
%   the function to do what you want, you might consider adapting these to your data.
%   
%   WindowLength    : Length of the windows (in seconds) for which correlation is computed; ideally
%                     short enough to reasonably capture periods of global artifacts or intermittent 
%                     sensor dropouts, but not shorter (for statistical reasons). Default: 5.
% 
%   MaxBrokenTime : Maximum time (either in seconds or as fraction of the recording) during which a 
%                   retained channel may be broken. Reasonable range: 0.1 (very aggressive) to 0.6
%                   (very lax). The default is 0.4.
%
%   NumSamples : Number of RANSAC samples. This is the number of samples to generate in the random
%                sampling consensus process. The larger this value, the more robust but also slower 
%                the processing will be. Default: 50.
%
%   SubsetSize : Subset size. This is the size of the channel subsets to use for robust reconstruction, 
%                as a fraction of the total number of channels. Default: 0.25.
%
% Out:
%   Signal : data set with bad channels removed
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2014-05-12

% Copyright (C) Christian Kothe, SCCN, 2014, christian@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License as published by the Free Software Foundation; either version 2 of the
% License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
% even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with this program; if not,
% write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
% USA

if ~exist('corr_threshold','var') || isempty(corr_threshold), corr_threshold = 0.8; end
if ~exist('noise_threshold','var') || isempty(noise_threshold), noise_threshold = 4; end
if ~exist('window_len','var') || isempty(window_len), window_len = 5; end
if ~exist('max_broken_time','var') || isempty(max_broken_time), max_broken_time = 0.4; end
if ~exist('num_samples','var') || isempty(num_samples), num_samples = 50; end
if ~exist('subset_size','var') || isempty(subset_size), subset_size = 0.25; end
if ~exist('reset_rng','var') || isempty(reset_rng), reset_rng = true; end

subset_size = round(subset_size*size(signal.data,1)); 

% flag channels
if max_broken_time > 0 && max_broken_time < 1  %#ok<*NODEF>
    max_broken_time = size(signal.data,2)*max_broken_time;
else
    max_broken_time = signal.srate*max_broken_time;
end

signal.data = double(signal.data);
[C,S] = size(signal.data);
window_len = window_len*signal.srate;
wnd = 0:window_len-1;
offsets = 1:window_len:S-window_len;
W = length(offsets);

fprintf('Scanning for bad channels...\n');

if signal.srate > 100
    % remove signal content above 50Hz
    B = design_fir(100,[2*[0 45 50]/signal.srate 1],[1 1 0 0]);
    parfor c = 1:signal.nbchan
        X(:,c) = filtfilt_fast(B,1,signal.data(c,:)'); 
    end

    % determine z-scored level of EM noise-to-signal ratio for each channel
    noisiness = mad(signal.data'-X)./mad(X,1);
    znoise = (noisiness - median(noisiness)) ./ (mad(noisiness,1)*1.4826);    

    % trim channels based on that
    noise_mask = znoise > noise_threshold;
else
    X = signal.data';
    noise_mask = false(C,1)'; % transpose added. Otherwise gives an error below at removed_channels = removed_channels | noise_mask';  (by Ozgur Balkan)
end

if ~(isfield(signal.chanlocs,'X') && isfield(signal.chanlocs,'Y') && isfield(signal.chanlocs,'Z') && all([length([signal.chanlocs.X]),length([signal.chanlocs.Y]),length([signal.chanlocs.Z])] > length(signal.chanlocs)*0.5))
    error('clean_channels:bad_chanlocs','To use this function most of your channels should have X,Y,Z location measurements.'); end

% get the matrix of all channel locations [3xN]
[x,y,z] = deal({signal.chanlocs.X},{signal.chanlocs.Y},{signal.chanlocs.Z});
usable_channels = find(~cellfun('isempty',x) & ~cellfun('isempty',y) & ~cellfun('isempty',z));
locs = [cell2mat(x(usable_channels));cell2mat(y(usable_channels));cell2mat(z(usable_channels))];
X = X(:,usable_channels);
  
% caculate all-channel reconstruction matrices from random channel subsets   
if reset_rng
    rng('default')
end
if exist('OCTAVE_VERSION', 'builtin') == 0
    P = hlp_microcache('cleanchans',@calc_projector,locs,num_samples,subset_size);
else
    P = calc_projector(locs,num_samples,subset_size);
end
corrs = zeros(length(usable_channels),W);
        
% calculate each channel's correlation to its RANSAC reconstruction for each window
timePassedList = zeros(W,1);
for o=1:W
    tic; % makoto
    XX = X(offsets(o)+wnd,:);
    YY = sort(reshape(XX*P,length(wnd),length(usable_channels),num_samples),3);
    YY = YY(:,:,round(end/2));
	corrs(:,o) = sum(XX.*YY)./(sqrt(sum(XX.^2)).*sqrt(sum(YY.^2)));
    timePassedList(o) = toc; % makoto
    medianTimePassed = median(timePassedList(1:o));
    fprintf('clean_channel: %3.0d/%d blocks, %.1f minutes remaining.\n', o, W, medianTimePassed*(W-o)/60); % makoto
end
        
flagged = corrs < corr_threshold;
        
% mark all channels for removal which have more flagged samples than the maximum number of
% ignored samples
removed_channels = false(C,1);
removed_channels(usable_channels) = sum(flagged,2)*window_len > max_broken_time;
removed_channels = removed_channels | noise_mask';

% apply removal
if mean(removed_channels) > 0.75
    error('clean_channels:bad_chanlocs','More than 75%% of your channels were removed -- this is probably caused by incorrect channel location measurements (e.g., wrong cap design).');
elseif any(removed_channels)
    try
        signal = pop_select(signal,'nochannel',find(removed_channels));
    catch e
        if ~exist('pop_select','file')
            disp('Apparently you do not have EEGLAB''s pop_select() on the path.');
        else
            disp('Could not select channels using EEGLAB''s pop_select(); details: ');
            hlp_handleerror(e,1);
        end
        fprintf('Removing %i channels and dropping signal meta-data.\n',nnz(removed_channels));
        if length(signal.chanlocs) == size(signal.data,1)
            signal.chanlocs = signal.chanlocs(~removed_channels); end
        signal.data = signal.data(~removed_channels,:);
        signal.nbchan = size(signal.data,1);
        [signal.icawinv,signal.icasphere,signal.icaweights,signal.icaact,signal.stats,signal.specdata,signal.specicaact] = deal([]);
    end
    if isfield(signal.etc,'clean_channel_mask')
        signal.etc.clean_channel_mask(signal.etc.clean_channel_mask) = ~removed_channels;
    else
        signal.etc.clean_channel_mask = ~removed_channels;
    end
end



% calculate a bag of reconstruction matrices from random channel subsets
function P = calc_projector(locs,num_samples,subset_size)
%stream = RandStream('mt19937ar','Seed',435656);
rand_samples = {};
for k=num_samples:-1:1
    tmp = zeros(size(locs,2));
    subset = randsample(1:size(locs,2),subset_size);
%    subset = randsample(1:size(locs,2),subset_size,stream);
    tmp(subset,:) = real(sphericalSplineInterpolate(locs(:,subset),locs))';
    rand_samples{k} = tmp;
end
P = horzcat(rand_samples{:});


function Y = randsample(X,num)
Y = [];
while length(Y)<num
    pick = round(1 + (length(X)-1).*rand());
    Y(end+1) = X(pick);
    X(pick) = [];
end
% 
% function Y = randsample(X,num,stream)
% Y = [];
% while length(Y)<num
%     pick = round(1 + (length(X)-1).*rand(stream));
%     Y(end+1) = X(pick);
%     X(pick) = [];
% end

function Y = mad(X,flag) %#ok<INUSD>
Y = median(abs(bsxfun(@minus,X,median(X))));

function B = design_fir(N,F,A,nfft,W)
% Design an FIR filter using the frequency-sampling method.
% The frequency response is interpolated cubically between the specified
% frequency points.function B = design_fir(N,F,A,nfft,W)
%   N : order of the filter
%   F : vector of frequencies at which amplitudes shall be defined
%       (starts with 0 and goes up to 1; try to avoid too 
%        sharp transitions)
%   A : vector of amplitudes, one value per specified frequency
%   nFFT : optionally number of FFT bins to use
%   W : optionally the window function to use (default: Hamming)
%   B : designed filter kernel
%
%  Christian Kothe, Swartz Center for Computational Neuroscience, UCSD, 2013

if nargin < 4 || isempty(nfft)
    nfft = max(512,2^ceil(log(N)/log(2))); end
if nargin < 5
    W = 0.54 - 0.46*cos(2*pi*(0:N)/N); end

% calculate interpolated frequency response
F = interp1(round(F*nfft),A,(0:nfft),'pchip');

% set phase & transform into time domain
F = F .* exp(-(0.5*N)*sqrt(-1)*pi*(0:nfft)./nfft);
B = real(ifft([F conj(F(end-1:-1:2))]));

% apply window to kernel
B = B(1:N+1).*W(:)';


function X = filtfilt_fast(varargin)
% Uses FFT convolution (needs fftfilt). The function is faster than filter when approx.
% length(B)>256 and size(X,Dim)>1024, otherwise slower (due size-testing overhead).
% Christian Kothe, Swartz Center for Computational Neuroscience, UCSD, 2010-07-14

if nargin == 3
    [B, A, X] = deal(varargin{:});
elseif nargin == 4
    [N, F, M, X] = deal(varargin{:});
    B = design_fir(N,F,sqrt(M)); A = 1; % note: we use the sqrt() because we run forward and backward
else
    help filtfilt_fast;
    return;
end

if A == 1
    was_single = strcmp(class(X),'single');
    w = length(B); t = size(X,1);    
    % extrapolate
    X = double([bsxfun(@minus,2*X(1,:),X(1+mod(((w+1):-1:2)-1,t),:)); X; bsxfun(@minus,2*X(t,:),X(1+mod(((t-1):-1:(t-w))-1,t),:))]);
    % filter, reverse
    X = filter_fast(B,A,X); X = X(length(X):-1:1,:);
    % filter, reverse
    X = filter_fast(B,A,X); X = X(length(X):-1:1,:);
    % remove extrapolated pieces
    X([1:w t+w+(1:w)],:) = [];
    if was_single
        X = single(X); end    
else    
    % fall back to filtfilt for the IIR case
    X = filtfilt(B,A,X);
end

function [X,Zf] = filter_fast(B,A,X,Zi,dim)
% Like filter(), but faster when both the filter and the signal are long.
% [Y,Zf] = filter_fast(B,A,X,Zi,Dim)
%
% Uses FFT convolution. The function is faster than filter when approx. length(B)>256 and
% size(X,Dim)>1024, otherwise slower (due size-testing overhead).
%
% See also:
%   filter, fftfilt
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2010-07-09
%
%                           contains fftfilt.m from Octave:
%                           Copyright (C) 1996-1997 John W. Eaton


% Copyright (C) Christian Kothe, SCCN, 2010, ckothe@ucsd.edu
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License as published by the Free Software Foundation; either version 2 of the
% License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
% even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with this program; if not,
% write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
% USA


if nargin <= 4
    dim = find(size(X)~=1,1); end
if nargin <= 3
    Zi = []; end

lenx = size(X,dim);
lenb = length(B);
if lenx == 0
    % empty X
    Zf = Zi;
elseif lenb < 256 || lenx<1024 || lenx <= lenb || lenx*lenb < 4000000 || ~isequal(A,1)
    % use the regular filter
    if nargout > 1
        [X,Zf] = filter(B,A,X,Zi,dim);
    else
        X = filter(B,A,X,Zi,dim);
    end
else
    was_single = strcmp(class(X),'single');
    % fftfilt can be used
    if isempty(Zi)
        % no initial conditions to take care of
        if nargout < 2
            % and no final ones
            X = unflip(oct_fftfilt(B,flip(double(X),dim)),dim);
        else
            % final conditions needed
            X = flip(X,dim);
            [dummy,Zf] = filter(B,1,X(end-length(B)+1:end,:),Zi,1); %#ok<ASGLU>
            X = oct_fftfilt(B,double(X));
            X = unflip(X,dim);
        end
    else
        % initial conditions available
        X = flip(X,dim);
        % get a Zi-informed piece
        tmp = filter(B,1,X(1:length(B),:),Zi,1);
        if nargout > 1
            % also need final conditions
            [dummy,Zf] = filter(B,1,X(end-length(B)+1:end,:),Zi,1); %#ok<ASGLU>
        end
        X = oct_fftfilt(B,double(X));
        % incorporate the piece
        X(1:length(B),:) = tmp;
        X = unflip(X,dim);
    end
    if was_single
        X = single(X); end
end

function X = flip(X,dim)
if dim ~= 1
    order = 1:ndims(X);
    order = order([dim 1]);
    X = permute(X,order);
end

function X = unflip(X,dim)
if dim ~= 1
    order = 1:ndims(X);
    order = order([dim 1]);
    X = ipermute(X,order);
end


function y = oct_fftfilt(b, x, N)
% Copyright (C) 1997 John W. Eaton
% If N is not specified explicitly, we do not use the overlap-add
% method at all because loops are really slow.  Otherwise, we only
% ensure that the number of points in the FFT is the smallest power
% of two larger than N and length(b). This could result in length
% one blocks, but if the user knows better ...

transpose = (size(x,1) == 1);
if transpose, x = x.'; end
[r_x,c_x] = size(x);
[r_b,c_b] = size(b);
if min([r_b, c_b]) ~= 1
    error('octave:fftfilt','fftfilt: b should be a vector'); 
end

l_b = r_b*c_b;
b = reshape(b,l_b,1);

if nargin == 2
    % Use FFT with the smallest power of 2 which is >= length (x) +
    % length (b) - 1 as number of points ...
    N = 2^(ceil(log(r_x+l_b-1)/log(2)));
    B = fft(b,N);
    y = ifft(fft(x,N).*B(:,ones(1,c_x)));
else
    % Use overlap-add method ...
    if ~isscalar(N)
        error ('octave:fftfilt','fftfilt: N has to be a scalar'); end
    N = 2^(ceil(log(max([N,l_b]))/log(2)));
    L = N - l_b + 1;
    B = fft(b, N);
    B = B(:,ones(c_x,1));
    R = ceil(r_x / L);
    y = zeros(r_x, c_x);
    for r = 1:R
        lo = (r - 1) * L + 1;
        hi = min(r * L, r_x);
        tmp = zeros(N, c_x);
        tmp(1:(hi-lo+1),:) = x(lo:hi,:);
        tmp = ifft(fft(tmp).*B);
        hi = min(lo+N-1, r_x);
        y(lo:hi,:) = y(lo:hi,:) + tmp(1:(hi-lo+1),:);
    end
end

y = y(1:r_x,:);
if transpose
    y = y.'; end

% Final cleanups: if both x and b are real respectively integer, y
% should also be
if isreal(b) && isreal(x)
    y = real(y); end
if ~any(b - round(b))
    idx = ~any(x - round(x));
    y(:,idx) = round(y(:,idx));
end

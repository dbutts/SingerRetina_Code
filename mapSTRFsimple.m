function [STRF,dt] = mapSTRFsimple( stim, spkts, frame_res, STIMRES, T, ss, silent )
%
% Usage: [STRF,dt] = mapSTRFsimple( stim, spkts, frame_res, <STIMRES>, <T>, <spatial_subset> )
%
% General STRF mapping with temporal resolution 'dt' with stim of the
% form stim(t,pixels)
%
% frame_res: time resolution of STRF in (STIMRES/frame_res) units
%            if not entered, default is 16 (~1 ms resolution)
%
% spatial_subset of form [xmin xmax ymin ymax]
%
% Modified from mapSTRF.m
% Created 28 September 2005 DAB
% Modified 11 October 2005 DAB

[STIMNUM NPIX] = size(stim);
if STIMNUM == 1
  stim = stim';  %'
  [STIMNUM NPIX] = size(stim);
end
if nargin < 7
	silent = 0;
end

% Minimum length of STRF (it will round to an integer amount of stimulus frames)
if NPIX > 1
  STRF_TIME = 400;  % ms
else
  STRF_TIME = 250;  % ms
end
if (nargin > 4)
  STRF_TIME = T*1000;
end

if nargin < 3
  frame_res = 8;  % Default ~1 ms resolution 
end

%%%%%%%%%%%%%%%%%%% Set constants and setup %%%%%%%%%%%%%%%%%%
if nargin < 4
  STIMRES = 8.33911056439 /1000; % 7.76363983 (for mseq)
else
  if STIMRES == 0
    STIMRES = 8.33911056439 /1000; % 7.76363983 (for mseq)
  end
end

%disp(sprintf( '%f %d %f', STIMRES, length(spkts), spkts(1) ))

% Set internal constants
NFRAMES = ceil(STRF_TIME/1000/STIMRES);
dt = STIMRES/frame_res;
NT = ceil(frame_res*NFRAMES);
KERNTIME = NFRAMES*STIMRES;
LPIX = sqrt(NPIX);

Nreps = length(find(diff(spkts) < 0));
if ~isempty(spkts) && (spkts(end) >= 0)
  Nreps = Nreps+1;
end

% Determine maximum stimulus time
T = STIMNUM*STIMRES;

spkts = spkts(find((spkts >= 0) & (spkts < T))); 

%FiringRate = length(spkts)/(STIMNUM*STIMRES);
FiringRate = length(spkts)/(size(stim,1)*STIMRES);

if ~silent
	disp(sprintf( '\nProcessing %d spikes.', length(spkts) ))
	disp(sprintf( '--------' ));
	disp(sprintf( '%5.2f Hz', FiringRate/Nreps ))
	disp(sprintf( '--------' ));
end

%%%%%%%%%%%%%%%%%%%% Calculate STRF %%%%%%%%%%%%%%%%%%%%%
if (nargin < 6) || isempty(ss)
  sprange = 1:NPIX;
else
  sprange = [];
  Lx = ss(2)-ss(1)+1;
  disp(sprintf( 'Output of %d x %d sub-spatial RF', Lx, ss(4)-ss(3)+1 ))
  for i=ss(3):ss(4)
    sprange(length(sprange)+(1:Lx)) = LPIX*(i-1)+(ss(1):ss(2));
  end
end

NPInc = length(sprange);

STRF = zeros(NT,NPInc);
if isempty(spkts)
	return
end

% Get relevant spikes
ts = spkts(find((spkts >= KERNTIME) & (spkts <= (STIMNUM*STIMRES))));
Nspks_used = length(ts);

for spkN = 1:length(ts)

  frameN = floor(ts(spkN)/STIMRES);                  % 0 is first frame
  subframeN = floor((ts(spkN)-frameN*STIMRES)/dt);   % 0 is first subframe

  indx = 1+frameN + floor((subframeN-(0:NT-1))/frame_res);  % stim index upsampled to dt

  STRF = STRF + stim(indx,sprange);
end

% Normalize to make 'average stimulus to evoke a spike'
if Nspks_used > 0
  STRF = STRF/Nspks_used;
end

%%%%%%%%%%%%%%%%%%%% Little analysis at the end %%%%%%%%%%%%%%%%%%
if NPIX > 1
  [maxrf maxpixel] = max(max(STRF));
  [minrf minpixel] = min(min(STRF));
  [maxrf maxt] = max(STRF(:,maxpixel));
  [minrf mint] = min(STRF(:,minpixel));
	if ~silent
		if maxt < mint  % then likely an ON cell
			disp(sprintf( '\nON cell: (max = %0.3f), lat = %0.2f ms, Pix #%d\n', maxrf, maxt*dt*1000, maxpixel ))
		else % likely off cell
			disp(sprintf( '\nOFF cell: (max = %0.3f), lat = %0.2f ms, Pix #%d\n', -minrf, mint*dt*1000, minpixel ))
		end
	end
else
  [maxrf maxt] = max(STRF);
  [minrf mint] = min(STRF);
	if ~silent
		if maxt < mint  % then likely an ON cell
			disp(sprintf( '\nON cell: (max = %0.3f), lat = %0.2f ms', maxrf, maxt*dt*1000 ))
		else % likely off cell
			disp(sprintf( '\nOFF cell: (max = %0.3f), lat = %0.2f ms', -minrf, mint*dt*1000 ))
		end
	end
end

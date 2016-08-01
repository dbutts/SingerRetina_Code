function newspks = clip_repeats( spks, repeat_list, trange )
%
% Usage: newspks = clip_repeats( spks, repeat_list, <trange> )
%
% Created long time ago DAB

spks = spks(:)';

if nargin < 3
  trange(1) = min(spks);
  trange(2) = max(spks);
  if trange(1) < 0
    trange(1) = 0;
  end
end

% Count number of trials (and positions in train)
[N locs] = repeat_detect(spks);

%disp(sprintf( 'Initially %d repeats.', N ))

rrange = repeat_list(repeat_list <= N);

newspks = [];

for i = 1:length(rrange)
  a = locs(rrange(i),1):locs(rrange(i),2);
  newspks(end+(1:length(a)+1)) = [spks(a) -1];
end

newspks = newspks(newspks < trange(2));
if trange(1) > 0
  newspks = newspks((newspks < 0) | (newspks > trange(1)));
end

function [fitstruct,Robs,modstim,Uindx,XVindx] = RETfit_LN( stim_params, Ndata, cc, blocks, trimstim, skipReg )
%
% Usage: [fitstruct,Robs,modstim,Uindx,XVindx] = RETfit_LN( stim_params, Ndata, cc, blocks_to_include, trimstim, skipReg )
%
% Use STA to trim stimulus around RF center (10 on each side), and then performs separable-LN fit to
% generate fitstructure with 

if nargin < 4
	if ~isempty(Ndata.blocks)
		if iscell(Ndata.blocks)
			blocks = Ndata.blocks{cc};
		else
			blocks = Ndata.blocks;
		end
	else
		blocks = [];
	end
end

[stim,spks] = format_noise_data( Ndata, cc, blocks );
fitstruct.cellname = Ndata.cellname{cc};

NT = size(stim,1);

L = 10;   % 21 on a side from max value

if (nargin < 5) || isempty(trimstim)
	trimstim = 1;
end
if nargin < 6
	skipReg = 0;
end

mod_params = stim_params;
modstim = stim;

frac = stim_params.up_fac;
NX = stim_params.dims(3);
NY = stim_params.dims(2);

% Format spike data
Robs = NIM.Spks2Robs( spks, stim_params.dt, NT*frac );

% Spike-triggered average
sta = sNIM.spike_triggered_average( Robs, stim, stim_params );

% Look at STA in time
if trimstim
	[~,PP] = max(max(abs(sta)));
	y = mod((PP-1),NY)+1;
	x = floor((PP-1)/NY)+1;
	x1 = max([x-L 1]);
	y1 = max([y-L 1]);
	x2 = x1+2*L;
	if x2 > 40
		x2 = 40; x1 = 40-2*L;
	end
	y2 = y1+2*L;
  if y2 > 30
		y2 = 30; y1 = 30-2*L;
	end
	indx = [];
	for x = x1:x2
		indx = [indx (x-1)*NY+(y1:y2)];
	end
	modstim = stim(:,indx);
	sta = sta(:,indx);
	mod_params.dims(2:3) = [1 1]*(2*L+1);
end

if size(sta,1) > stim_params.dims(1)
	sta = sta(1:mod_params.dims(1),:);
else
	sta(end:mod_params.dims(1),:) = 0;
end

% Initialize linear (separable) model
fitstruct.LLs = [];
fitstruct.SPindx = indx;

sLN = sNIM( sta, 1, mod_params);
sLN = sLN.set_reg_params( 'd2x', 10, 'd2t', 10 );
[Uindx,XVindx] = sLN.generate_XVfolds( size(stim,1), 5, 3 );

% Fit Separable LN model
sLN = sLN.fit_TSalt( Robs, modstim, Uindx, 'silent', 0 );

% optimize regularization 
if ~skipReg
	sLN = sLN.reg_pathT( Robs, modstim, Uindx, XVindx );
	sLN = sLN.reg_pathSP( Robs, modstim, Uindx, XVindx );
end
sLN = sLN.correct_spatial_signs();

fitstruct.Uindx = Uindx;
fitstruct.XVindx = XVindx;
[LL0,~,~,fitprops] = sLN.eval_model( Robs, modstim, XVindx );
fitstruct.LLs = LL0-fitprops.nullLL;
fitstruct.LN = sLN;

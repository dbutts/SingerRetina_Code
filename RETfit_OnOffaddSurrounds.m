function newstruct = RETfit_OnOffaddSurrounds( fitstruct, Ndata, cc, skipReg, manual_class )
%
% Usage: newstruct = RETfit_OnOffaddSurrounds( fitstruct, Ndata, cc, <skipReg>, <manual_class> )
%
% Default is to skip regularization, since fit will be refined from here

if (nargin < 4) || isempty(skipReg)
	skipReg = 1;
end
if nargin < 5
	if isfield(fitstruct,'class')
		manual_class = fitstruct.class;
	else
		manual_class = [];
	end
end

frac = fitstruct.LN.stim_params.up_fac;
[stim,spks,dt] = format_noise_data( Ndata, cc );
modstim = stim(:,fitstruct.SPindx);
Robs = NIM.Spks2Robs( spks, fitstruct.LN.stim_params.dt, size(stim,1)*frac );
%Robs = histc( spks, 0:(dt/frac):(size(stim,1)*dt) ); 
%Robs = Robs(1:size(stim,1)*frac);

[newstruct,autoclass] = RETfit_classifyOnOff( fitstruct );
%newstruct.class = autoclass;
if ~isempty(manual_class)
	autoclass = manual_class;
	disp( 'Using manual classification.' )
end
if autoclass == 0
	disp('Class = 0. Not fitting.')
	return 
end	
	
if autoclass == 2
	disp( 'Identified On-Off.')
	fit0 = fitstruct.OnOff;
	clone = [1 2];
else
	if autoclass == 1
		disp( 'Identified On.')
	else
		disp( 'Identified Off.')
	end
	fit0 = fitstruct.LN;
	clone = 1;
end
for nn = 1:length(clone)
	fit0.subunits(clone(nn)).NLtype = 'rectlin';
	fit0.subunits(clone(end)+nn) = fit0.subunits(nn);
	fit0.subunits(clone(end)+nn).NLtype = 'rectlin';
    fit0.subunits(clone(end)+nn).kt = fit0.subunits(nn).timeshift_kts(3); %temporal kernel delay 
	fit0.subunits(clone(end)+nn).weight = -1; %suppressive unit 
end

fit1 = fit0.fit_TSalt( Robs, modstim, fitstruct.Uindx, 'fit_offsets', 1 );
fprintf( '  Add sup: %f -> %f\n',fit0.fit_props.LL, fit1.fit_props.LL ); 

if ~skipReg
	fit1 = fit1.reg_pathT( Robs, modstim, fitstruct.Uindx, fitstruct.XVindx, 'subs', length(clone)*2:-1:1 );
	fit1 = fit1.reg_pathSP( Robs, modstim, fitstruct.Uindx, fitstruct.XVindx, 'subs', length(clone)*2:-1:1 );
	fit1 = fit1.reg_pathT( Robs, modstim, fitstruct.Uindx, fitstruct.XVindx );
end

[LLnew,~,~,fitprops] = fit1.eval_model( Robs, modstim, fitstruct.XVindx );
newstruct.LLs(end+1) = LLnew-fitprops.nullLL;
newstruct.SupFit = fit1;
fprintf( '           XV: %f -> %f\n', max(fitstruct.LLs), newstruct.LLs(end) )


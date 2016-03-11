function fitstruct = RETfit_LN2OnOff( LNstruct, Ndata, cc, skipReg )
%
% Usage: fitstruct = RETfit_LN2OnOff( LNstruct, Ndata, cc, <skipReg> )
%
% Default is to skip regularization, since fit will be refined from here

if nargin < 4
	skipReg = 1;
end
fitstruct = LNstruct;

% format data
frac = LNstruct.LN.stim_params.up_fac;
[stim,spks,dt] = format_noise_data( Ndata, cc );
modstim = stim(:,LNstruct.SPindx);
Robs = histc( spks, 0:(dt/frac):(size(stim,1)*dt) ); 
Robs = Robs(1:size(stim,1)*frac);

% Build ON-OFF
OOfit = LNstruct.LN;
OOfit.subunits(1).NLtype = 'rectlin';
OOfit.subunits(2) = OOfit.subunits(1);
Rect3fit = OOfit;
OOfit.subunits(2).kt = -OOfit.subunits(1).kt;
OOfit = OOfit.fit_TSalt( Robs, modstim, LNstruct.Uindx, 'silent', 1 );
OOfit = OOfit.correct_spatial_signs();
fprintf( 'Cell %d:\tLN -> OnOff: %f -> %f\nRegularizing-T...\n', cc, LNstruct.LN.fit_props.LL, OOfit.fit_props.LL ); 
OOfit = OOfit.reg_pathT( Robs, modstim, LNstruct.Uindx, LNstruct.XVindx, 'subs', [2 1], 'silent', 1 );

if ~skipReg
	OOfit = OOfit.reg_pathSP( Robs, modstim, LNstruct.Uindx, LNstruct.XVindx, 'subs', [2 1] );
	OOfit = OOfit.reg_pathT( Robs, modstim, LNstruct.Uindx, LNstruct.XVindx );
end

fitstruct.OnOff = OOfit;
fitstruct.LLs(2) = OOfit.eval_model( Robs, modstim, LNstruct.XVindx );
fprintf( '   XV: \t%f -> %f\n', fitstruct.LLs(1), fitstruct.LLs(2) )

return;

%% Rect-3 stuff is only useful if not ON-OFF
Rect3fit.subunits(2).weight = -1;
Rect3fit.subunits(3) = Rect3fit.subunits(2);
Rect3fit.subunits(3).kt = -Rect3fit.subunits(3).kt;
Rect3fit.subunits(2).kt = LNstruct.LN.subunits(1).timeshift_kts(3);
Rect3fit = Rect3fit.fit_TSalt( Robs, modstim, LNstruct.Uindx, 'silent', 1 );
Rect3fit = Rect3fit.correct_spatial_signs();
fprintf( '        \tOnOff -> Rect3: %f -> %f\n',OOfit.fit_props.LL, Rect3fit.fit_props.LL ); 
if ~skipReg
	Rect3fit = Rect3fit.reg_pathT( Robs, modstim, LNstruct.Uindx, LNstruct.XVindx, 'subs', [2 1] );
	Rect3fit = Rect3fit.reg_pathSP( Robs, modstim, LNstruct.Uindx, LNstruct.XVindx, 'subs', [2 1] );
end
fitstruct.Rect3 = Rect3fit;
fitstruct.LLs(3) = Rect3fit.eval_model( Robs, modstim, LNstruct.XVindx );
fprintf( '   XV: \t\t\t-> %f\n', fitstruct.LLs(3) )


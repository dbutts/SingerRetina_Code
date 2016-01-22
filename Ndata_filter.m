function NdataMod = Ndata_filter( Ndata, NblockMod )
%
% Usage: NdataMod = Ndata_filter( Ndata, NblockMod )
%
% filter Ndata to only include Cells/Blocks given by NblockMod, which is generated by Ndata_valid_blocks

for ctype = 1:length(Ndata)
	if ~isempty(Ndata{ctype})
		Ncells = 0;
		for cc = 1:length(Ndata{ctype}.cellname)
			blocks = NblockMod{ctype}{cc};
			if blocks(1) > 0
				Ncells = Ncells + 1;
				NdataMod{ctype}.cellname{Ncells} = Ndata{ctype}.cellname{cc};
				NdataMod{ctype}.Quals(cc) = Ndata{ctype}.Quals(cc);
				NdataMod{ctype}.stiminfo{Ncells}{length(blocks)} = [];
				NdataMod{ctype}.spks{Ncells}{length(blocks)} = [];
				NdataMod{ctype}.blocks{cc} = blocks;
				for nn = 1:length(blocks)
					NdataMod{ctype}.stiminfo{Ncells}{nn} = Ndata{ctype}.stiminfo{cc}{blocks(nn)};
					NdataMod{ctype}.spks{Ncells}{nn} = Ndata{ctype}.spks{cc}{blocks(nn)};
				end
			end
		end
		NdataMod{ctype}.dt = Ndata{ctype}.dt;
	end
end

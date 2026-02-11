@markovjunior 'b' begin
	# Build a scaffold of center areas and boundaries
    @do_n 5 begin
        @rule b => G
    end
	@block repeat begin
		@do_n 1 begin
			@rule Gbb => GYY
		end
		@do_all begin
			@sequential
			@rule YYR => BBR
			@rule BBb => RRb
			@rule YYb => BYY
			@rule BYY => BBR
		end
	end

	# Build tiles around that scaffold.
	@block repeat begin
		@do_n 300 begin
			@rule bb => gg
			@infer begin
				@path 5 bB => bBM => R
			end
		end
		@do_all begin
			@rule gb => gM
		end
	end

	# Clean up the scaffold and other bits.
	@do_all begin
		@sequential
		# Remove R and expand G.
		@rule R => B
		@rule BGB => GGG

		# Relax B into either M or g.
		@rule gBg => ggg
		@rule gBM => ggM
		# Any Blue left is now on the edges.
		@rule gB => gg
		@rule MB => MM

		# Shore up the M walls.
		@rule b => M
		@rule MgM => MMM
	end

	# Use Green areas as a seed for low points.
	@do_all begin
		@sequential
		@rule Gg => GG
	end
	# Fill in medium-elevation areas with Brown (for mountains).
	@do_all begin
		@sequential
		@rule GMg => GMN
		@rule GMMg => GMMN
		@rule GMMMg => GMMMN
		@rule Ng => NN
	end
	# Fill in high-elevation areas with white (for peaks).
	@do_all begin
		@rule g => w
	end

	# Prune Magenta borders.
	@do_all begin
		# Don't mess with low areas as much; Magenta can be used for water
		# @rule GMG => GGG
		# @rule GMMG => GGGG
		@rule NMN => NNN
		@rule NMMN => NNNN
		@rule wMw => www
		@rule wMMw => wwww
	end
	# Interpret Magenta differently based on where it is.
	@do_all begin
		@sequential
		# Between low and medium areas, it's a cliff face.
		@rule GMN => GEN
		@rule GMMN => GEEN
		@rule EM => EE
	end
	# Around low areas, it's a forest or water.
	@block repeat begin
		@do_n 1 begin
			@rule GMG => GLG
			@rule GMMG => GLLG
			@rule GMMMG => GLLLG
			@rule GMG => GBG
			@rule GMMG => GBBG
			@rule GMMMG => GBBBG
		end
		@do_all begin
			@rule LM => LL
			@rule BM => BB
		end
	end
	# Finally, any leftover Magenta takes on a neighbor's value.
	@do_all begin
		@sequential
		@rule GM => GG
		@rule NM => NN
		@rule wM => ww
		@rule M => L
	end
end





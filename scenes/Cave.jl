@markovjunior 'I' begin
	# Place a seed
    @do_n 1 begin
        @rule I => b
    end

	# Randomly carve outward, and
	#    randomly carve outward with partial (grey) rock
    @do_n 5000 begin
        @rule bI => gb
        @rule gI => bg
        @rule gb => bb
		@rule gb => gg
    end

	# Clean up the grey rock a bit
	@do_n 100 begin
		@rule bgb => bbb
	end
	@do_all begin
		@rule Ig => II
		@rule Ig => gg
		@rule Ig => bb
	end

	# Finalize colors.
	@do_all begin
		@rule g => S
	end

	# Add walls around the edges
	@do_all begin
		@rule bI => bS
	end
end
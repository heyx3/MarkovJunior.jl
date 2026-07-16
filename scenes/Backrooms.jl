@markovjunior 'I' begin
	# Generate a regular old maze.
	@rewrite 1 I=>Y
	@rewrite YII => YEY
	@rewrite Y => E

	# Mangle it
	@rewrite I=>E %(0.3:0.5)

	# Clean up the walls into something coherent.
	@rewrite begin
		PRIORITIZE(earliest)
		IEI => III  %0.25 \[ x ]
		IEI => III        \[ y ]
	end
	@rewrite [
		E E E
		E I E
		E E E
	  ] => [
		_ _ _
		_ E _
		_ _ _
	]	 \[ (x, y)[ (+x, +y) ] ]

	# Ensure areas are well-connected.
	#   1) Pick a home region
	@rewrite 1 E=>M
	@rewrite ME => MM
	#   2) Keep connecting it to other regions
	@sequence repeat begin
		# The first op of a 'repeat' sequence determines whether it should quit.
		# We want to quit when the whole map is connected.
		@rewrite 1 E => T
		@rewrite 1 T => E

		# Grab any low-hanging fruit.
		@rewrite 1 MIIIIIE => MMMMMME
		@rewrite 1 MIIIIE => MMMMME
		@rewrite 1 MIIIE => MMMME
		@rewrite 1 MIIE => MMME
		@rewrite 1 MIE => MME
		@rewrite ME => MM

		# Now try a loop-erased random walk along the boundary of the home region.
		@rewrite 1 MI => MR
		@rewrite (area/10) begin
			PRIORITIZE(earliest)
			 RE => LE  # Make it to a destination!
			RIE => LEE # Make it to a destination!

			OGB => IIR # Finish the loop-erasure
			OGG => IIO # Continue the loop-erasure

			RII => GGR # Continue the walk
			RIG => OIB # Nowhere else to go; give up and start erasing
		end
		# If there's an L on the grid, then we got a new path and need to carve it out.
		# Otherwise this search was a dud and needs to be deleted.
		@rewrite L[ROBG] => LL
		@rewrite [ROBG] => I
		@rewrite M[LE] => MM
	end
	# Give up on everything else.
	@rewrite E=>I
	@rewrite M=>E
end
@markovjunior begin
	# Generate a regular old maze.
	@rewrite 1 b=>Y
	@rewrite Ybb => YEY
	@fill 'E' +Y

	# Mangle it.
	@rewrite b=>E %(0.3:0.5)

	# Clean up the walls into something coherent.
	@rewrite begin
		PRIORITIZE(earliest)
		bEb => bbb  %0.25 \[ x ]
		bEb => bbb        \[ y ]
	end
	@rewrite [
		E E E
		E b E
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
		@rewrite 1 MbbbbbE => MMMMMME
		@rewrite 1 MbbbbE => MMMMME
		@rewrite 1 MbbbE => MMMME
		@rewrite 1 MbbE => MMME
		@rewrite 1 MbE => MME
		@rewrite ME => MM

		# Now try a loop-erased random walk along the boundary of the home region.
		@rewrite 1 Mb => MR
		@rewrite (area/10) begin
			PRIORITIZE(earliest)
			 RE => LE  # Make it to a destination!
			RbE => LEE # Make it to a destination!

			OGB => bbR # Finish the loop-erasure
			OGG => bbO # Continue the loop-erasure

			Rbb => GGR # Continue the walk
			RbG => ObB # Nowhere else to go; give up and start erasing
		end
		# If there's an L on the grid, then we got a new path and need to carve it out.
		# Otherwise this search was a dud and needs to be deleted.
		@rewrite L[ROBG] => LL
		@rewrite [ROBG] => b
		@rewrite M[LE] => MM
	end
	# Finalize colors.
	@fill 'E' +M

	# Add lights.
	#=
	@rewrite 1 E=>Y
	@sequence repeat begin
		@rewrite area/40 E=>L
		@rewrite area/180 E=>Y
	end field(LY->E)
	@rewrite L=>E
	=#
end
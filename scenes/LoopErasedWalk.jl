@markovjunior begin
	# This is a building block for more complex algorithms.
	# It builds random depth-first paths through the grid,
	#    and when it gets trapped it rewinds a mostly-random amount.
	# This version runs->rewinds->runs forever, but you can easily have it stop
	#    by replacing the Red head with a special "done" color
	#    when you get to where you're going.

	@rewrite 1 b=>R
	@rewrite begin
		PRIORITIZE(earliest)
		OGB => bbR # Finish the loop-erasure
		OGG => bbO # Continue the loop-erasure
		Rbb => GGR # Continue the walk
		RbG => ObB # Nowhere else to go; give up and start erasing
	end
end
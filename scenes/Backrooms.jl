@markovjunior begin
	# Generate a maze
	@rewrite 1 b=>Y
	@rewrite Ybb => YEY
	@fill 'E' +Y

	# Mangle the maze
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
	] \[ (x, y)[ (+x, +y) ] ]

	# Pick some areas, flood-fill them, and remove everything else.
	@rewrite (area/8000) E=>G
	@rewrite GE => GG
	@fill 'b' +E

	# Try to connect the two areas if they aren't already,
	#   using another flood-fill to detect that.
	# If they're too far apart they'll just stay separate
	#   (we need to implement a Path Op to 100% ensure connection).
	#  1) flood-fill one area:
	@rewrite 1 G=>R
	@rewrite [RE]G => _E
	@rewrite R => E
	#  2) try to connect any paths between that and the other area:
	@rewrite (area/2000) begin
		PRIORITIZE(earliest)
		GbbbbE => EEEEEE
		GbbbE => EEEEE
		GbbE => EEEE
		GbE => EEE
	end
	#  3) Clean up colors
	@fill 'E' +G
end
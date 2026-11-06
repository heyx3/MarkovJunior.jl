@markovjunior 2 begin
	# This one is a proof-of-concept, for how
	#   real games can work with MarkovJunior.jl.
	# It generates a detailed level,
	#   emits an event for the user code to copy out that level,
	#   then builds a minimap of that same level.
	
	# Generate some terrain.
	@rewrite 5 b => w
	@rewrite (area/50) wb => ww
	@rewrite wb => wN
	@rewrite Nb => NN %0.45
	@rewrite Nb => NG
	@rewrite ((area/25):(area/3)) Gb => GG
	@rewrite Gb => GL
	@rewrite ((area/25):(area/5)) Lb => LL
	@rewrite b => B
	@rewrite (area/15) L => E field(B -> L)
	@rewrite (area/20) B => I field(-E)

	# Place some human settlements.
	@rewrite 3 G => M
	@rewrite 3 N => O

	# Export the full map.
	@event full_map

	# Now generate a mini-map, preserving the important landmarks.
	@downscale 2... begin
		# If any of the important landmarks are in a downscaled chunk,
		#   use their color.
		M <= has(M)
		O <= has(O)
		# Otherwise, use the most popular color in each chunk.
		common
	end
end
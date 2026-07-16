@markovjunior begin
	@pragma GuiMaterial B metal 0.2 (0.1, 0.1, 0.8)

	# Place three seeds.
	@rewrite 1 b => R
	@rewrite 1 b => G
	@rewrite 1 b => B
	# Prevent any touching seeds.
	@rewrite begin
		R[GB]b => Rb[2]
		GBb => GbB
		bGB => GbB
	end

	# Grow the seeds.
	# Wherever they touch, recoil using Magenta.
	@rewrite begin
		# Start a recoil:
		R[GB] => MM  *100
		GB => MM     *100
		# Grow/close a recoil:
		M[RGB] => MM  *10
		M[RGB] => MO  *8
		M      => O   *2.5
		# Otherwise focus on growing the seeds:
		[RGB]b => [1][1]  /10
	end

	# Empty out the recoil space.
	@rewrite O => b
	# Smooth out the blobs.
	@rewrite b[RGB]b => bbb
end
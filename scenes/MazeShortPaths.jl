@markovjunior 'b' begin
	# This is a building block for more complex algorithms.
	# It generates a maze by building breadth-first paths.

    @rewrite 1 b => w
	@rewrite bbw => wGw
	@rewrite G => w
end
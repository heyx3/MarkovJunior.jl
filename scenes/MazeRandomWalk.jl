@markovjunior 'b' begin
	# This is a building block for more complex algorithms.
	# It generates a maze by building depth-first paths.
	# As soon as it gets stuck somewhere it rewinds until it finds the next opening.

    @rewrite 1 b => R
    @rewrite begin
		PRIORITIZE(earliest)
        Rbb => GGR
        GGR => Rww
    end
	@rewrite R => w
end
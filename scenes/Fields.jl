@markovjunior begin
	# A basic field will bias rewrite actions
	#   to be closer to source cells;
	#   in this case white ones.
    @rewrite 1 b=>w
    @rewrite b=>g field(w)
	@rewrite w=>g

	# Fields can be flipped, and randomized.
	@rewrite 1 g=>R
	@rewrite g=>G field(-R, random)
	@rewrite R=>G

	# The randomization can be finely controlled.
	@rewrite 1 G=>B
	@rewrite G=>w field(B, randomness=0.2)
	@rewrite B=>w

	# The fields can be constrained to specific path cells.
	#   1) generate a maze of w' paths and 'b' walls
    @fill 'b'
	@fill 'w' uv(center=0.5, size=0)
	@rewrite wbb=>wgw
	@fill 'w' +g
	#   2) place seeds and fill the spaces outward from those seeds.
	@rewrite (length/12) w=>P
    @rewrite (area/6) w=>O field(P -> w)
    @rewrite (area/6) w=>L field(O -> w)
	@rewrite          w=>N field(L -> w)
    @fill 'O' +P
    # The above can be sort-of accomplished with a flood-fill rewrite op,
    #   but that op couldn't ensure the coverage is even --
    #   it might go really far down one way and not much down another!

	# There are other features too, which are showcased
    #    with other, more interesting sample scenes!
end
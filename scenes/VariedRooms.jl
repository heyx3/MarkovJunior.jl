@markovjunior 2 begin
	# Generates some rooms, then uses the 'within_regions' Op
	#   to fill each room very differently.

	# Build a funky room plan.
	@rewrite begin
		[
		  b b b
		  b b b
		  b b w
		 ] => [
		  _ _ _
		  _ w w
		  _ w w
		 ]  %0.25  /40
		[
		  b b b b
		  b b b b
		  b w w b
		 ] => [
		  _ _ _ _
		  _ w w _
		  _ w w _
		 ]  %0.45  /10

		[
		  b b b b
		  b b b b
		  b b b b
		  b b b b
		 ] => [
		  _ _ _ _
		  _ w w _
		  _ w w _
		  _ _ _ _
		]  %0.15   /100

		[
		  b b b b b b
		  b b b b b b
		  b w w w b b
		 ] => [
		  _ _ _ _  _   _
		  _ w w w  _   _
		  _ w w w {wb} _
		]  %0.25   /5

		[
			_ w _
			w b w
			_ w _
		 ] => [
			_ _ _
			_ w _
			_ _ _
		]  /1000 # Generally only do it after all the others

		wbw => www %0.2
		wbbw => wwww %0.1
	end

	# Upscale so we can fill in detail.
	@upscale 4...

	# Detail each room individually.
	#   1) Assign a random color to each, indicating room-type.
	@sequence repeat begin
		@rewrite 1 w => {RGB}
		@rewrite w[RGB] => [2]_
	end
	#   2) Mark the boundaries using Yellow.
	@rewrite b[RGB] => bY
	#   3) Add the detail!
	@within_regions begin
		# Red rooms: pillars around the wall, carpet on the inside
		RY => @sequence begin
			# Stagger pillars around the edges.
			@rewrite YYY => YIY
			
			# Add carpet through the middle.
			@rewrite begin
				[
				  [Rg] [Rg] [Rg]
				   R    R    R
				  [Rg] [Rg] [Rg]
				 ] => [
				  _ _ _
				  g g g
				  _ _ _
				]
				gRg => _g_
			end
			@rewrite gI => R_
			@rewrite [
				I _
			    _ g
			  ] => [
				_ _
				_ R
			]
		end

		# Green rooms: grid of pillars
		GY => @sequence begin
			@rewrite [
				Y  Y    _
				Y  G   [YG]
				_ [YG]  _
			] => [
				_ _ _
				_ I _
				_ _ _
			]
			@rewrite [
				[YG] [YG] [YG] [YG] [YG]
				[YG]  I    G    G   [YG]
				[YG] [YG] [YG] [YG] [YG]
  			] => [
				_ _ _ _ _
				_ I G I _
				_ _ _ _ _
			]
			#@rewrite IGG => ILI
			@rewrite L => G
		end

		# Blue rooms: desks, books, and people using them
		BY => @sequence begin
			@rewrite YB => N_  %0.1
			@rewrite NN => NY

			@rewrite [
				B B B B B
				B B B B B
				B B B B B
			  ] => [
				g  g   {gE}  g   g
				g  N    N    N   g
				g {gE}  g   {gE} g
			] %0.2
			@rewrite NB => _E  %0.3
		end
	end
				
	# Clean up.
	@rewrite [RGB] => w
	@rewrite Y => w
end
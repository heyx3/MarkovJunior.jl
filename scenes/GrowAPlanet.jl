@markovjunior 2 'b' begin
	# Add fractal detail through a repeated sequence.
	@sequence 3 begin

		# Place some seeds.
		@rewrite (area/200) b=>g

		# Double the grid size along each axis.
		@upscale 2...

		# Grow the seeds into irregular blobby shapes.
		@rewrite (area/10) begin
			gb => gg
			[ g g
			  g b
			] => [
			  _ _
			  _ g
			]   *3  # stronger weight makes the shapes more rectangular
		end
	end

	# Let the gray drip downward.
	@rewrite (area/100) gbbb => gggb \[ +y ]

	# Add color.
	@rewrite bbg => bbL \[ +y ]
	@rewrite (area/400) L => G
	@rewrite (area/200) LG => GG
	@fill 'S' +b
	@fill 'N' +g
	@rewrite GN => GL
	@rewrite LN => LL %0.25
	
	# Add clouds.
	@rewrite 5 S => w
	@rewrite (area/10) begin
		wS => ww
		[ w w
		  w S
		] => [
		  _ _
		  _ w
		]  *5
	end
	# Emboss the clouds to separate them
	#    from the sky behind and ground in front.
	@rewrite Sw => Sg \[ -x, -y ]
	@rewrite [NGL]w => _b \[ +x, +y ]
	@rewrite bw => bg
end
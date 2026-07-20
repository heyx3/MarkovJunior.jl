@markovjunior 2 'I' begin
	# Place cave seeds.
	@fill 'R' uv(min=0.1, size=0)
	@fill 'R' uv(min=0.5, size=0)
	@fill 'R' uv(min=0.9, size=0)

	# Grow the cave seeds and leave some veins behind.
    @rewrite (area*4.5) begin
        R[bI] => bR
		IRb => G__
    end

	# Clean up the veins and reduce their frequency.
	@rewrite begin
		PRIORITIZE(earliest)
		R => b
		GG => SS  %0.1
		G => b
    end

	# Turn the veins into real minerals.
	@sequence repeat begin
		# Mark a vein as either Gold or Nitra.
		@rewrite 1 begin
			S => R  *3
			S => Y
		end
		# Flesh out that vein.
		@rewrite begin
			[RY]S => _[1]
		end
	end
end
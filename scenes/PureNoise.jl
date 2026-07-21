@markovjunior begin
	# A quick test of masks.
    @rewrite begin
		b    => w %0.5   # Masked to happen in ~half the grid
		[bw] => g %0.167 # Masked to only happen in ~1/3 of the places where the first rule applies
	end
end
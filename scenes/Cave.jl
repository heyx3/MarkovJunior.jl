@markovjunior 'I' begin
	# Place a seed.
	@fill 'b' uv(center=0.5, size=0)

	# Randomly carve outward, and
	#    randomly carve outward with partial (Slate) rock.
    @rewrite area*2.5 begin
        bI => Sb
		SI => bS
		Sb => bb
		Sb => SS
	end

	# Add a dripping grey shadow under the Slate rock.
	@rewrite begin
		Sb => gb \[+y]
		# Drip, and optionally forbid further drip using a Green pixel.
		gbbb => ggg{Gb} %0.4 \[+y]
	end
	@rewrite G => b
	@rewrite bg => bb  \[ +y ] # Remove drip that has no source

	# Add walls around the edges.
	@rewrite [bg]I => _S
end
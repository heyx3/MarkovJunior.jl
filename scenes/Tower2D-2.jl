@markovjunior 2 'T' begin
    # Draw guidelines that the tower edges must stay within,
    #   and initial drawing "heads" at the inside-bottom of either side of the tower.
    @fill 'R' uv(min=(0.2, 0.23), max=(0.4, 1))
    @fill 'G' uv(min=(0.4, 1), size=0)
    @fill 'R' uv(min=(0.6, 0.23), max=(0.8, 1))
    @fill 'B' uv(min=(0.6, 1), size=0)

    # Each drawing head has three colors --
    #   'drawing upward' (GB), 'drawing outward' (LI), and 'drawing inward' (YS).
    # We start upward.
    @rewrite begin
        # Normally just draw forwards.
        [GB]R => b[1]  \[ -y ]
        [IY]R => b[1]  \[ +x ]
        R[LS] => [2]b  \[ +x ]

        # Occasionally, indicate a "turn signal" instead.
        # Only turn in valid directions.
        R[GB] => R[LS]  /8  \[ +x ]
        [GB]R => [YI]R  /8  \[ +x ]
        [LI] => [GB]    /4
        [YS] => [GB]    /4
    end
    # Clean up the guidelines.
    @fill 'T' +R
    # Connect the two sides at the top.
	@rewrite [LYIS] => [GGBB]
	@rewrite begin
		PRIORITIZE(earliest)
		GB => bb
		G[Tb] => bG \[ +x ]
		B[Tb] => bB \[ -x ]
	end

    # Fill the inside of the tower.
    @fill 'g' uv(center=0.5, size=0)
    @rewrite gT => gg
	# Trim hanging walls.
	@rewrite TbT => TTT
	
	# Place potential window markers.
	#  1) Add a white border to prevent windows from touching the outside
	@rewrite bg => bw
	@rewrite [
		  b w
		  w g
  	] => [
  		_ _
  		_ w
  	] \[ (x, y)[ (x, y) ] ]
	#  2) Mark each inside corner Red
	@rewrite [
  		w g
  		w w
  	] => [
  		_ R
  		_ _
  	]  \[ (x, y)[ (x, y) ] ]
	#  3) draw the inner corners outward to find 1D connections,
	#       roughly segmenting the building surface.
	#     Each innner corner has two possible directions to search,
	#        so transition it to Blue and then Orange after each one.
	@sequence repeat begin
		# Start a search and transition its source cell
		@rewrite 1 begin
			# Look for horizontal connections first,
			#   to avoid vertical ones from reaching all the way to bottom
			# PRIORITIZE(earliest)
			w[RB]g => w[BO]Y  *10  \[ x ]
			w[RB]g => w[BO]Y  \[ y ]
		end

		@rewrite 1 [BO]Y[RGBO] => _G_   # Handle short connections
		@rewrite 1 [BO]Yg => __M  # Mark the search direction
		@rewrite begin
			PRIORITIZE(earliest)

			# Normally advance forward.
			YMg => YYM
			
			# If we hit an interesting block, cap off our search
 		   #   with Olive to indicate success.
			YM[RBOGw] => LL_
			# If we hit a problematic block/botto of the building,
			#   remove the Magenta to indicate failure.
			YM => _Y

			# Spread success/failure across the line.
			LY => LL
		end
		# Lock in success as Green, and failure as background grey
		@rewrite begin
			L => G
			Y => g
		end
	end
	# Convert the inner corners to just be part of the Green lines.
	@rewrite [RBO] => G
	#  4) Place windows in each area within the green lines.
	#     Exclude the bottom part of the tower.
	@fill 'M' uv(center=(0.5, 1), size=0)
	@rewrite Mg => MM
	#  5) Clean up colors.
	@fill 'Y' +g
	@fill 'g' +MGw

	# Lastly, add antennae to the top of the building.
	@rewrite [
  	  T T T
	    b b b
      ] => [
		T R T
		b b b
	  ]  %0.2 \[ (x, y)[ (+x, +y) ] ]
	@rewrite RT => RR %0.1 \[ -y ]
end
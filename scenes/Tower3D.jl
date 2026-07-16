@markovjunior 3 begin
	@pragma fast_fills

    # Finalized colors:
    #   L: grass/leaves
	@pragma GuiMaterial L dielectric  0.75 (0, 0.5, 0.2)
    #   G: bushes
	@pragma GuiMaterial G dielectric  0.475 (0, 1, 0)
    #   N: bark/floors/tables
	@pragma GuiMaterial N dielectric  0.75 (0.5, 0.2, 0)
    #   E: cubicle walls
	@pragma GuiMaterial E dielectric  1.0 (1, 0.9, 0.8)
    #   Y: lamps
	@pragma GuiMaterial Y light_source   (10, 10, 6)
    #   g: furniture/ceiling
	@pragma GuiMaterial g dielectric  1.0  0.5
	#   I: windows/guidelines
	@pragma GuiMaterial I glass  0.4  1.0 (0, 0.2, 0.5)

	# Define the lot
	@fill 'L' uv(min=0, max=(1, 1, 0))
	@fill 'g' uv(min=(0.27, 0.27, 0), max=(0.73, 0.73, 0))

	# Define the bounding box of the building.
	@fill 'I' uv(min=(0.27, 0.27, 0), max=(0.27, 0.73, 1))
	@fill 'I' uv(min=(0.27, 0.27, 0), max=(0.73, 0.27, 1))
	@fill 'I' uv(min=(0.27, 0.73, 0), max=(0.73, 0.73, 1))
	@fill 'I' uv(min=(0.73, 0.27, 0), max=(0.73, 0.73, 1))

	# Plant bushes
	@rewrite [
		L L L L L
		L L L L L
		L L L L L
		L L L L L
		L L L L L ;;;
		b b b b b
		b b b b b
		b b b b b
		b b b b b
		b b b b b
	  ] => [
		L L L L L
		L L L L L
		L L L L L
		L L L L L
		L L L L L ;;;
		b b b b b
		b b b b b
		b b G b b
		b b b b b
		b b b b b
	  ]   %0.1  \[ x[ +x ], y[ +y ], z[ +z ] ]

	# Turn some bushes into trees
	@rewrite G => N  %0.2
	@rewrite LNb => LNR
	@rewrite LNRb => LNRR
	@rewrite NRb => NNR  %0.2
	@rewrite [RB]b => wb \[ +z ]
	@sequence repeat begin
		# Add one leaf layer to each tree
		@rewrite [
			b b b
			b w b
			b b b
		  ] => [
			L L L
			L Y L
			L L L
		  ]  \[ x[ +x ], y[ +y ] ]
		# Decide whether each tree should get another layer
		@rewrite begin
			Yb=>Mb  /10
			Yb=>Lw
		end
        @rewrite Y => L # Cap off trees that hit the ceiling
		# Round off the tips of trees that are finished
		@rewrite [
			L L L
			L M L
			L L L
		  ] => [
			b L b
			L L L
			b L b
		  ] \[ x[ +x ], y[ +y ] ]
	end
	@rewrite w=>L
	@rewrite R=>N

	# Build upwards
	@fill 'w' uv(center=(0.5, 0.5, 0), size=0)
	@sequence repeat begin
		# Start a new floor.
		# Each floor has 4 z-levels, with different logic on each
		@rewrite 1 wbbbb => RBMTw  \[ +z ]

        # The first level is a wood floor
        @rewrite Rb => RR \[ x, y ]
		# Randomly shrink the floorspace a bit in one direction
		#TODO: Implement once we have some meta-ops like @rewrite_line
		#=
		@rewrite 1 [
			I I I I
			R R R R
		  ] => [
			I I I I
			P O O P
		]    %0.01   \[ x[x,y],  y[x,y]  ]
		@rewrite 1:3 [
			P O O P
			R R R R
		  ] => [
			P O O P
			P O O P
		]   \[ x[x,y], y[x,y]  ]
	    @rewrite OPR => SOP  \[ x, y ]
		@rewrite [SO] => P
		@rewrite Pbbb => gIII \[ +z ]
		@rewrite Ib => II \[ +z ]
		=#
        @rewrite R=>N

        # The second level adds walls and furniture
        @rewrite Bb => BB \[ x, y ]
        # The third level distinguishes between full walls and cubicle walls, and adds lamps
        @rewrite Mb => MM \[ x, y ]
        # The third level will be ceiling and ceiling-lights.
        @rewrite Tb => TT \[ x, y ]

        #   * Use a scaled-up maze generator to make areas
        @rewrite IB => IO   # Rim of empty space around the room
        @rewrite 1 [  # Seed the maze at a corner
            I I I I
            I O O O
            I O B B
            I O B B
          ] => [
            I I I I
            I O O O
            I O O O
            I O O O
          ]   \[ x[x,y], y[x,y] ]
        @rewrite [
            O O   B B   B B
            O O   B B   B B
          ] => [
            O O   R R   O O
            O O   R R   O O
          ]   \[  x[x,y],  y[x,y]  ]
        #   * Open up the areas a bit
        @rewrite RBBR => ROOR %0.2 \[ x, y ]
        @rewrite R => O
        # Add furniture to each area
        @rewrite [
            B B B
            O O O
            O O O
          ] => [
            _ _ _
            _ N _
            _ _ _
        ]  %0.15 \[ (x, y)[ (+x, +y), (+y, +x), (x, -y), (-y, x) ] ]
        @rewrite [
            _    B    B    B
            [OB] O    O    B
            O    O    O    B
            [OB] [OB] [OB] _
          ] => [
            _ _ _ _
            _ g g _
            _ _ g _
            _ _ _ _
        ]  \[ (x, y)[ (+x, +y), (+y, +x), (x, -y), (-y, x) ] ]
        # Finalize colors.
        @rewrite [OR] => b
        @rewrite B => E

        # Fill in the third level
        #   * Add lamps or empty space above each desk
        @rewrite begin
            NNM => NNY  /4  \[ +z ]
            NNM => NNb       \[ +z ]
        end
        #   * Add empty space above empty space
        @rewrite NbM => __b   \[ +z ]
        #TODO: pick cubicle walls vs full walls
        @rewrite EM => EE \[ +z ]
        @rewrite M => b

        # Fill in the fourth level with ceiling tiles and lights
        @rewrite ET => EE \[ +z ]  # Extend walls upward to help light-placement logic.
        @rewrite ETTTE => EgYgE  %0.4  \[ x, y ]
        @rewrite ETTTTE => EgYggE  %0.2  \[ x, y ]
		@rewrite ETTE => EYTE  %0.1 \[ x, y ]
        @rewrite T => g
	end
	# Put a ceiling at the top
	@rewrite w => T
	@rewrite Tb => TT \[ x, y ]
	@rewrite T => g

	# Replace the guideline with outer walls
	@rewrite I => g

	# Add external doors.
	@rewrite 1:2 [
		L g g
		b g b
		b g b
	  ] => [
		_ _ _
		_ b _
		_ b _
	]  \[  x[x,y],  y[+z]  ]
	#TODO: Add an elevator covering one strip along the outside

	# Place windows around each level
	@rewrite [
		N g b
		b g b
		b g b
	  ] => [
		N g b
		b {gI} b
		b I b
	]  %(0.1:0.3)  \[ x[x, y], y[ +z ] ]
end
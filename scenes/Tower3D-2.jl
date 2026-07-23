@markovjunior 3 'T' begin
	@pragma fast_fills
	@pragma GuiMaterial b dielectric 0.85  0.01
	@pragma GuiMaterial g dielectric 0.15  0.5
	@pragma GuiMaterial w dielectric 0.85  1
	@pragma GuiMaterial T empty

    # Draw a road next to the building's future location
	@fill 'G' uv(min=(0, 0, 0), max=(1, 0.07, 0))
	@fill 'b' uv(min=(0, 0.07, 0), max=(1, 0.35, 0))
	@fill 'g' uv(center=(0.5, 0.19, 0), size=(1, 0, 0))
	@fill 'w' uv(center=(0, 0.19, 0), size=0)
	@rewrite wgg => wbw
	# Add guard-rails
	@rewrite [
		b
		[TG] ;;;
		T
		T
	  ] => [
		g
		_ ;;;
		g
		_
	]  \[ (x, y, z)[ (+x, y, +z) ] ]

	# Draw a parking lot around the building.
	@fill 'b' uv(min=0, max=(1, 1, 0)) +T
	@fill 'R' uv(min=(1, 1, 0), size=0)
	@rewrite 2 Rbb=>ROB
	@rewrite OBbbbb => OBbbOB
	@rewrite [
		b b
		B O
	  ] => [
		b w
		b w
	]  \[ x[x, y], y[x, y] ]
	@rewrite begin
		wwb => RRR *10  \[ x, y ]
		wwb => bbb      \[ x, y ]
	end
	@rewrite [ROB] => w

	@rewrite [
		b b b b b b
		b b b b b b
		b b b b b b
		b b b b b b
		b b b b b b
		b b b b b b
      ] => [
		_ _ _ _ _ _
		_ _ _ _ _ _
		_ _ R R _ _
		_ _ R R _ _
		_ _ _ _ _ _
		_ _ _ _ _ _
	]  %0.5 \[ x[+x], y[+y] ]
	@sequence (length/3) begin
		@rewrite begin
			# Most of the time, draw another layer upward:
			[
				R R
				R R ;;;
				T T
				T T
			  ] => [
				R R
				R R ;;;
				G G
				G G
			]   *10  \[ (x, y, z)[ (+x, +y, +z) ] ]

			# Sometimes, stop drawing early.
			#TODO: To feasibly support this we need a 'directional' bias; disabled for now
			[
				R R
				R R ;;;
				T T
				T T
			  ] => [
				R R
				R R ;;;
				B B
				B B
			] *0  \[ (x, y, z)[ (+x, +y, +z) ] ]
		end
		@rewrite [
			G G
			G G
		  ] => [
			R R
			R R
		] \[ (x, y)[ (+x, +y) ] ]
	end
	@rewrite [ R;R ;; R;R ;;; T;T ;; T;T ] => [ g;g ;; g;g ;;; _;_ ;; _;_ ] \[ (x, y)[ (+x, +y) ] ]
	@rewrite [ R;R ;; R;R ;;; g;g ;; g;g ] => [ g;g ;; g;g ;;; _;_ ;; _;_ ] \[ (x, y)[ (+x, +y) ] ]

	# Draw a connecting layer between all the "stilts"
	
end




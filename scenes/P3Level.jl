@markovjunior 'b' 3 begin
	#=
		Materials:

			b: solid floor (in-between playable floors)
			I: solid wall
			R: room open area
			G: room doorway
			B: room vertical connection
			Y: room hallway
			w: unprocessed solid floor
			L: in-process solid floor
            gMNOPETS: short-term use
	=#
	@pragma GuiMaterial I glass  0.4  1.0 (0, 0.2, 0.5)
	@pragma GuiMaterial w empty

	# Mark some floors as reachable and some as solid space.
	@fill 'w' uv(min=0, size=(1, 1, 0))
	@rewrite wbb => wbw \[ z ]

    @sequence repeat begin
		# Flood one floor with L
	    @rewrite 1 w => L
		@rewrite Lw => LL \[ x, y ]
		
		# Generate rooms.
		@sequence repeat begin
			# Place some basic areas.
			@rewrite begin
				[ L _ _ _ L
				  _ L L L _
                  _ L L L _
				  _ L L L _
				  L _ _ _ L
				] => [
				  _ _ _ _ _
				  _ N N N _
				  _ N N N _
				  _ N N N _
				  _ _ _ _ _
				]  %0.03   \[ (x, y)[ (+x, +y) ] ]
				NLN => _G_  %(0.2:0.2)  \[ x, y ]
			end

			# Refine the output of the previous step.
			@rewrite begin
				# Remove lone square rooms.
				[ L L L L L
				  L N N N L
				  L N N N L
				  L N N N L
				  L L L L L
				] => [
				  L L L L L
				  L L L L L
				  L L L L L
				  L L L L L
				  L L L L L
			    ] \[ (x, y)[ (+x, +y) ] ]
			end
		end

		# Carve hallways between disparate rooms.
		@rewrite 1 N => O
		@rewrite begin
			NO => OO \[ x, y ]
			NGO => O__ \[ x, y ]
		end
		@sequence repeat begin
			# Loop-erased walk to find a link between O and R
			@rewrite 1 NL => _P
			@rewrite (area/10) begin
				PRIORITIZE(earliest)
				 PO => T_     \[ x, y ] # Make it to a destination!
				PLO => MT_    \[ x, y ] # Make it to a destination!

				SME => PMM    \[ x, y ] # Finish the loop-erasure
				SMM => LLS    \[ x, y ] # Continue the loop-erasure

				PLL => MMP    \[ x, y ] # Continue the walk
				PLM => SLE    \[ x, y ] # Nowhere else to go; give up and start erasing
			end
			# If there's a T on the grid, then we got a new path and need to carve it out.
			# Otherwise this search was a dud and needs to be deleted.
			@rewrite TM => TT   \[ x, y ]
			@rewrite T => Y     \[ x, y ]
			@rewrite begin
				NO => O_   \[ x, y ]
				NGO => O__ \[ x, y ]
				YN => _O   \[ x, y ]
			end
			@rewrite [PSEM] => L     \[ x, y ]
		end

		# Finish up.
		@rewrite [LON] => [IRR]
	end

	# Find vertical connections between floors; mark them 'B'.
	@rewrite (length*0.7) RbR => RBR %0.01 \[ z ]

	@event Minimap

	# Stretch the floor-plan out into a voxel scene.
	@upscale (5, 5, 5)
end
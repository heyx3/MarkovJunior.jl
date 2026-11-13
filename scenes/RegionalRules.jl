@markovjunior begin
    # A neat algorithm that demonstrates how to run sequences within a sub-region of the grid.

	# Grow part-way out from the center of the level.
	@fill 'w' uv(center=0.5, size=texel)
	@rewrite (area/1.5) wb => ww

	# Draw some streaks back inward.
	@rewrite begin
		PRIORITIZE(earliest)

		Rw => RR %0.85 \[ -x ]
		Ow => OO %0.85 \[ +x ]
		Mw => MM %0.85 \[ -y ]
		Pw => PP %0.85 \[ +y ]

		bw => _R %0.25  \[ -x ]
		bw => _O %0.25  \[ +x ]
		bw => _M %0.25  \[ -y ]
		bw => _P %0.25  \[ +y ]
	end
	@rewrite [ROMP] => b

	# Plop a maze in the center.
	@within_box uv(center=0.5, size=0.4) begin
		@rewrite w => G
		@rewrite 1 G=>g
		@rewrite gGG => gRg
		@rewrite g => R
	end

	# If the maze has walls on the bottom row, put it on a pedestal.
	@rewrite Rw => _I \[ +y ]
	@sequence repeat begin
		# Randomly draw downward from the current bottom row.
		@rewrite Iw => IT  %0.4 \[ +y ]
		# Remove the "bottom row" status from that row and give it to the next.
		@rewrite I => B
		@rewrite T => I

		# Find the convex hull of those draws.
		@rewrite begin
			PRIORITIZE(earliest)
			# Put out N feelers rightward;
			#    if we hit another I then rewind those N's into more I's.
			
			NI => II \[ +x ]
			Nw => NN \[ +x ]

			Iw => IN \[ +x ]
		end
		@rewrite N => w
	end
	# It's quite likely the pedestal didn't reach the ground.
	@rewrite Iw => BI \[ +y ]
	@rewrite I => B

	# If the maze has walls on the top row,
	#    build upward streaks to the ceiling.
	@rewrite Rw => RO %0.2 \[ -y ]
	@rewrite ROw => RRO \[ -y ]

	# Finalize colors.
	@rewrite [RBTO] => b
	@rewrite G => w
end
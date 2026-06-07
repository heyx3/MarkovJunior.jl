@markovjunior 2 'T' begin
    @fill 'b' uv(min=(0.35, 1), max=(0.65, 1))

	# Build the skeleton of the tower.
	@sequence repeat begin
        # Draw upwards, occasionally placing markers for side rooms.
        # Remember that this first step needs to always have at least one application
        #   or else the entire sequence stops repeating.
		@rewrite begin
			PRIORITIZE(earliest)
			bT => RT %0.1 \[ x ]
			bT => wg \[ -y ]
		end

        # Draw side rooms at the marker spots.
		@rewrite [
			T T T
			T T T
			T T R
		] => [
			G G G
			E T L
			E E E
		] \[ x[ x ] ]
        @rewrite R => w  # Remove markers that don't fit

        # Start the next stage of drawing upwards
		@rewrite g => b

        # Occasionally place platform markers along the top edge.
		@rewrite bT => YT %0.06 \[ x ]
	end
    # Clean up tower's main colors
	@rewrite [bgR] => w

	# Draw out the platforms
	@rewrite begin
		# Y = cursor
		# N = platform
		# R = platform (above support)
		# I = support

		PRIORITIZE(earliest)

		# Draw the platform sideways; prevent supports (R) from touching
        YTT => NTT  %0.2  \[ x ]  # Cut the platform off sometimes
        RTT => RNY        \[ x ]  # Ensure supports can't touch
        YTT => NRT  %0.4  \[ x ]  # Sometimes place a platform
                                  #  (remember all rules in a @rewrite share the same mask RNG,
                                  #    so 0.2 here would never happen as
                                  #    the previous 0.2 rule always supercedes it)
		YTT => NYT        \[ x ]  # Normal behavior

		# Draw down the supports
		TRT => TNI \[ +y ]
		IT => II \[ +y ]

		# Clean up
		Y => N
		R => I
	end
	@rewrite w => I

    # Draw out the side rooms
	@rewrite G => T
	@rewrite IIL => ITT
	@rewrite L => T

    # Place some windows
    @rewrite begin
        PRIORITIZE(earliest)

        [
            I I I I I
            I I I I I
            I I I I I
			I I I I I
			I I I I I
        ] => [
			I    I     I     I   I
            I    I     Y     I   I
            I    Y     Y     Y   I
            I    I     Y     I   I
			I    I     I     I   I
        ]  %0.02

        [
            I I I I I
            I I I I I
            I I I I I
			I I I I I
			I I I I I
        ] => [
			I    I     I     I   I
            I    I     I     I   I
            I    Y     Y     Y   I
            I    I     I     I   I
			I    I     I     I   I
        ]  %0.05

        [
            I I I
            I I I
            I I I
        ] => [
            _ _ _
            _ Y _
            _ _ _
        ]  %0.07
    end

	# Create a night sky.
	@fill 'b' uv(min=0, max=(1, 0.85)) +T %0.975
	@fill 'b' uv(min=(0, 0.85), max=1) +T
	@fill 'g' +T
end
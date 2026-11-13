@markovjunior 'b' 2 begin
    # Mark the min corner.
    @fill 'B' uv(min=0, size=texel)
    # Pick an "across-brick" axis (ideally vertical but let's make it axis-agnostic for fun).
    @rewrite 1 Bbb => BYY

    # Mark the rows where each brick line starts.
    @rewrite begin
        YYbbbbbbb   => GGGGGGRYY
        YYbbbbbbbbb => GGGGGGGGRYY
    end
    # Push the marker to the end of the across-brick axis.
    @rewrite begin
        PRIORITIZE(earliest)
        YYb => GYY
        YY => GG
    end
    # Draw out each row.
    @rewrite begin
        PRIORITIZE(earliest)
        Bbb => BYY
        Rbb => RYY
        YYb => IYY
        YY => II
    end

    # Now fully draw out each row, starting with the first.
    @sequence repeat begin
        @rewrite BII => BYY
        # As it's drawn out, insert column markers for the bricks underneath.
        @rewrite begin
            YYIIIIIIIIIII     => TTTTTTTTTTMYY
            YYIIIIIIIIIIIII   => TTTTTTTTTTTTMYY
            YYIIIIIIIIIIIIIII => TTTTTTTTTTTTTTMYY
        end
        @rewrite YYI => TYY
        @rewrite YY => TT

        # Draw downward to the next row.
        #   1. Start with the min edge, and mark the row below it to eventually go through this same process.
        @rewrite BGG => OYY
        @rewrite begin
            PRIORITIZE(earliest)
            YYG => OYY
            YYR => OOB
            YY => OO # At the bottom of the grid there is no other row to process
        end
        #   2. Draw the rest of the column markers down.
		@rewrite Mbb => TYY
		@rewrite YYb => PYY
		@rewrite YY => PP
        #   3. Fill in the bricks.
        @rewrite begin
            PRIORITIZE(earliest)
            YYb => LYY
            YY => LL
            Tbb => TYY
        end
    end

    # Finalize the colors!
    # Start with larger lines to speed up the process.
    @rewrite OOOO => TTTT
    @rewrite OOO => TTT
    @rewrite OO => TT
    @rewrite O => T
	@rewrite TTTTTTTT => wwwwwwww
	@rewrite TTTTTTT => wwwwwww
	@rewrite TTTTTT => wwwwww
	@rewrite TTTTT => wwwww
	@rewrite TTTT => wwww
	@rewrite TTT => www
	@rewrite TT => ww
	@rewrite T => w
	@rewrite PPPP => gggg
	@rewrite PPP => ggg
	@rewrite PP => gg
	@rewrite P => g
	@rewrite LLLLLLLLLL => RRRRRRRRRR
	@rewrite LLLLLLLL => RRRRRRRR
	@rewrite LLLLLL => RRRRRR
	@rewrite LLLL => RRRR
	@rewrite LL => RR
	@rewrite L => R
end
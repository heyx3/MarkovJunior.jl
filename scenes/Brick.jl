@markovjunior 'b' 2 begin
    # Mark the min corner.
    @draw_box 'B' min=0 size=0

    # Pick an "across-brick" axis (ideally vertical but whatever).
    @do_n 1 begin
        @rule Bbb => BYY
    end
    # Mark the rows where each brick line starts.
    @do_all begin
        @rule YYbbbbbbb => GGGGGGRYY
        @rule YYbbbbbbbbb => GGGGGGGGRYY
    end
    # Push the marker to the end of the across-brick axis.
    @do_all begin
        @sequential
        @rule YYb => GYY
        @rule YY => GG
    end
    # Draw out each row.
    @do_all begin
        @sequential
        @rule Bbb => BYY
        @rule Rbb => RYY
        @rule YYb => IYY
        @rule YY => II
    end

    # Now fully draw out each row, starting with the first.
    @block repeat begin
        @do_all begin
            @rule BII => BYY
        end
        # As it's drawn out, insert column markers for the bricks underneath.
        @do_all begin
            @rule YYIIIIIIIIIII => TTTTTTTTTTMYY
            @rule YYIIIIIIIIIIIII => TTTTTTTTTTTTMYY
            @rule YYIIIIIIIIIIIIIII => TTTTTTTTTTTTTTMYY
        end
        @do_all begin
            @sequential
            @rule YYI => TYY
            @rule YY => TT
        end

        # Draw downward to the next row.
        #   1. Start with the min edge, and mark the row below it to eventually go through this same process.
        @do_all begin
            @rule BGG => OYY
        end
        @do_all begin
            @sequential
            @rule YYG => OYY
            @rule YYR => OOB
            @rule YY => OO # At the bottom of the grid there is no other row to process
        end
        #   2. Draw the rest of the column markers down.
        @do_all begin
            @sequential
            @rule Mbb => TYY
            @rule YYb => PYY
            @rule YY => PP
        end
        #   3. Fill in the bricks.
        @do_all begin
            @sequential
            @rule YYb => LYY
            @rule YY => LL
            @rule Tbb => TYY
        end
    end

    # Finalize the colors!
    @do_all begin
        # Start with larger lines to speed up the process.
        @sequential

        @rule OOOO => TTTT
        @rule OOO => TTT
        @rule OO => TT
        @rule O => T

        @rule TTTTTTTT => wwwwwwww
        @rule TTTTTTT  => wwwwwww
        @rule TTTTTT   => wwwwww
        @rule TTTTT    => wwwww
        @rule TTTT     => wwww
        @rule TTT      => www
        @rule TT       => ww
        @rule T        => w

        @rule PPPP => gggg
        @rule PPP  => ggg
        @rule PP   => gg
        @rule P    => g

        @rule LLLLLLLLLL => RRRRRRRRRR
        @rule LLLLLLLLL => RRRRRRRRR
        @rule LLLLLLLL => RRRRRRRR
        @rule LLLLLLL => RRRRRRR
        @rule LLLLLL => RRRRRR
        @rule LLLLL => RRRRR
        @rule LLLL => RRRR
        @rule LLL => RRR
        @rule LL => RR
        @rule L => R
    end
end
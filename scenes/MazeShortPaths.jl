@markovjunior 'b' begin
    # Pick a source cell
    @do_n 1 begin
        @rule "b" => "w"
    end
    # Draw maze paths out from the cell
    @do_all begin
        @rule "bbw" => "wEw"
    end
    # Clean up
    @do_all begin
        @rule "E" => "w"
    end
end
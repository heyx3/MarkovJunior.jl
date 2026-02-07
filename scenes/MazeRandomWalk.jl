@markovjunior 'b' begin
    # Pick a source pixel
    @do_n 1 begin
        @rule "b" => "R"
    end
    # Do a backtracking random walk to carve out maze paths
    @do_all begin
        @sequential
        @rule "Rbb" => "GGR"
        @rule "GGR" => "Rww"
        @rule "R" => "w"
    end
end
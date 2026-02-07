@markovjunior 'I' begin
    @do_n 1 begin
        @rule I => R
    end
    @do_n 5000 begin
        @rule R_ => bR
    end
    @do_all begin
        @rule R => b
    end
end
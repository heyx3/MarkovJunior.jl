@markovjunior 'I' begin
    @do_n 1 begin
        @rule I => b
    end
    @do_n 500 begin
        @rule bI => gb
        @rule gI => bg
        @rule gb => bb
    end
end
@markovjunior 'b' begin
    @do_n 5 begin
        @rule b => G
    end
    @do_all begin
        @sequential
        @rule Gbb => GBB
        @rule Gb => GB
    end
    @do_all begin
		@sequential
		@rule BbB => BRB
        @rule BBb => BBB
    end
	@do_all begin
		@rule bb => YY
		@infer begin
			@path recompute 0 R => bB => G
		end
	end
end
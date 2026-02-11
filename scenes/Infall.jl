@markovjunior 'b' begin
    @do_n 5 begin
        @rule b => G
    end
	@block repeat begin
		@do_n 1 begin
			@rule Gbb => GYY
		end
		@do_all begin
			@sequential
			@rule YYR => BBR
			@rule BBb => RRb
			@rule YYb => BYY
			@rule BYY => BBR
		end
	end

	@do_all begin
		@rule bb => YY
		@infer begin
			@path 3 bB => bB => R
		end
	end
end
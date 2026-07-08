@markovjunior begin
	# Place seeds
    @rewrite (length/10) b=>w

	# Place tiles growing outward from that seed
	@rewrite wb => wg
	@sequence repeat begin
		#  1) grow the tile
		@rewrite ((area/8):(area/5.5)) gb=>gg field(w, random)
		#  2) draw the tile's border
		@rewrite bg => Mg
		#  3) initialize the next outer tile
		@rewrite Mb => Mg
	end
	# Break some of the borders to merge adjacent tiles
	@rewrite (length/15 : length/10) gMg=>ggg

	# Flood-fill from the seeds to create plains
	@rewrite w => G
	@rewrite Gg => GG

	# Add mountains/oceans next to plains,
	#    plains next to oceans,
	#    and peaks next to mountains
	@sequence repeat begin
		@rewrite 1 begin
			PRIORITIZE(earliest)
			GMg => GM{IN}
			IMg => IMG
			NMg => NMw
			wMg => wMw
		end
		@rewrite [NIGw]g => [1][1]
	end
	
	# Resolve the borders between biomes:
	#  * Between mountain and plains is beige
	@rewrite begin
		PRIORITIZE(latest)
		NMG => NEG
		EM => EE
	end
	#  * Between two plains is river/lake/forest
	@rewrite begin
		PRIORITIZE(latest)
		[BL]M => [1][1]
		[
			M   G
			G  [BL]
		] => [
		  [2,2]  G
			G    _
		]
		GMG => G{BL}G
	end
	#  * Between two oceans is more ocean
	@rewrite IMI => III
	#  * Between mountain and peak is more mountain and peak
	@rewrite NMw => N{Nw}w
	#  * Other borders randomly choose their neighbors
	@rewrite [GIN]M => [1][1]
	@rewrite [BLw]M => [1][1]
end





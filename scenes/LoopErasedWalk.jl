@markovjunior begin
	@rewrite 1 b=>R
	@rewrite begin
		PRIORITIZE(earliest)
		OGB => bbR # Finish the loop-erasure
		OGG => bbO # Continue the loop-erasure
		Rbb => GGR # Continue the walk
		RbG => ObB # Nowhere else to go; give up and start erasing
	end
end
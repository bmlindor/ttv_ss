using JLD2

function extract_data(filename)
	f = jldopen(String(filename), "r")
	return f
end
#unsure if this works within for loops
#might have used in 200000 step moon run on desktop

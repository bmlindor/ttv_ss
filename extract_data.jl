using JLD2

function extract_data(filename)
	f = jldopen(String(filename), "r")
end

function show_args(args)
@show args
end

show_args(ARGS)
x,y = ARGS[1],parse(Float64,ARGS[2])
#y = show_args(ARGS)[1]
#z = parse(Float64,show_args(ARGS)[2])
function f(x,y)
return y*2
end

println(typeof(x))
println(typeof(y))

z = f(x,y)

println(typeof(z))

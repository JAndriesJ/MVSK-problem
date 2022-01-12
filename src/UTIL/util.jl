module util
export get_ord

get_ord = x -> Float64(floor(Int, log10(abs(x))))

function run_tests()
    @assert get_ord(10^-12) == -12.0

end


end
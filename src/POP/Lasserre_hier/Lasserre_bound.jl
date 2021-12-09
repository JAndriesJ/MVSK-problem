module Lasserre_bound
include("SDPmodel.jl")
include("SDPoptimized.jl")
using .SDPmodel
using .SDPoptimized

export get_Lasserre_bound

function get_Lasserre_bound(Data_matrix,t)
    model = get_SDP_model(Data_matrix,t)
    return optimize_SDP(model)
end

end
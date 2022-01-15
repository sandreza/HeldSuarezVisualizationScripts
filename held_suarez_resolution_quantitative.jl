include("held_suarez_resolution.jl")
##
field_function_list = (
    (; comparison_field = "u", func = grab_state),
    (; comparison_field = "T", func = grab_state),
    (; comparison_field = "TT", func = eddy_variance),
    (; comparison_field = "uv", func = eddy_variance),
    (; comparison_field = "vT", func = eddy_variance),
    (; comparison_field = "uu", func = eddy_variance),
    (; comparison_field = "v", func = grab_state),
    (; comparison_field = "vv", func = eddy_variance),
    (; comparison_field = "wT", func = eddy_variance),
    (; comparison_field = "ρ", func = grab_state),
    (; comparison_field = "ρe", func = grab_state),
    (; comparison_field = "w", func = grab_state),
)

field_function_list = field_function_list[1:6]
error_matrix = zeros(5, length(field_function_list))
for (i, field_function) in enumerate(field_function_list)
    comparison_field = field_function.comparison_field
    grab = field_function.func

    jl_file = jldopen(filename[1], "r+")
    u_truth = grab(comparison_field, jl_file)
    u_truth isa Tuple ? u_truth = u_truth[1] : nothing

    jl_file = jldopen(filename[2], "r+")
    u_2 = grab(comparison_field, jl_file)
    u_2 isa Tuple ? u_2 = u_2[1] : nothing

    jl_file = jldopen(filename[3], "r+")
    u_3 = grab(comparison_field, jl_file)
    u_3 isa Tuple ? u_3 = u_3[1] : nothing

    jl_file = jldopen(filename[4], "r+")
    u_4 = grab(comparison_field, jl_file)
    u_4 isa Tuple ? u_4 = u_4[1] : nothing

    jl_file = jldopen(filename[5], "r+")
    u_5 = grab(comparison_field, jl_file)
    u_5 isa Tuple ? u_5 = u_5[1] : nothing

    jl_file = jldopen(filename[6], "r+")
    u_6 = grab(comparison_field, jl_file)
    u_6 isa Tuple ? u_6 = u_6[1] : nothing

    error_list = Float64[]
    er1 = abs.(u_truth - u_2) / maximum(abs.(u_truth))
    er2 = abs.(u_truth - u_3) / maximum(abs.(u_truth))
    er3 = abs.(u_truth - u_4) / maximum(abs.(u_truth))
    er4 = abs.(u_truth - u_5) / maximum(abs.(u_truth))
    er5 = abs.(u_truth - u_6) / maximum(abs.(u_truth))
    qn = 1.0 # quantile number
    push!(error_list, quantile(er1[:], qn))
    push!(error_list, quantile(er2[:], qn))
    push!(error_list, quantile(er3[:], qn))
    push!(error_list, quantile(er4[:], qn))
    push!(error_list, quantile(er5[:], qn))
    # push!(error_list, maximum(abs.(u_truth - u_7)))

    error_matrix[:, i] .= error_list

    minerror = title_names[argmin(error_list)+1] # since u_truth is 1
    println("The minimum error for " * comparison_field, " is ", minerror)
    println("All the errors are ", error_list)
    println("-------------------------")

end
# defined in utils.jl
##
include("utils.jl")
latex_format(error_matrix)

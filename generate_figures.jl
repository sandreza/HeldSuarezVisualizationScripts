figure_directory = "figures/"

# Discontinuous Galerkin Intro visualizations 
include("dg_projection.jl")
save(figure_directory * "DG_projections.png", fig)
include("dg_derivative.jl")
save(figure_directory * "DG_derivatives.png", fig)
include("dg_split_form_tendencies.jl")
save(figure_directory * "DG_split_form_tendencies.png", fig)

# Held-Suarez visualizations
include("held_suarez_statistics.jl")
save(figure_directory * "Held_Suarez_Statistics.png", fig)

include("held_suarez_resolution.jl")
save(figure_directory * "Held_Suarez_Resolutions.png", fig)

include("small_held_suarez_statistics.jl")
save(figure_directory * "Small_Held_Suarez_Statistics.png", fig)

include("held_suarez_snapshots.jl")
save(figure_directory * "surface_fields.png", fig)

include("held_suarez_3d_snapshots.jl")
save(figure_directory * "vertical_velocity.png", fig)

include("annulus_plot.jl")
save(figure_directory * "annulus_and_gridpoints.png", fig)

include("cubed_sphere_warp.jl")
save(figure_directory * "cubed_sphere.png", fig)
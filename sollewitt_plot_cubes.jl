using PlotlyJS
using NPZ

# --- CUBE GEOMETRY DEFINITION (1-based indexing) ---
vertices = [
    0 0 1; 0 1 1; 1 1 1; 1 0 1;  # Top face
    0 0 0; 0 1 0; 1 1 0; 1 0 0   # Bottom face
]

edges = [
    0 1; 1 2; 2 3; 3 0;  # Top edges
    4 0; 5 1; 6 2; 7 3;  # Vertical edges
    4 5; 5 6; 6 7; 7 4   # Bottom edges
] .+ 1 # Convert to 1-based indexing

# --- CONNECTIVITY FUNCTION ---
function is_connected(visibility_vector, edges_array)
    active_edge_indices = findall(==(1), visibility_vector)

    if isempty(active_edge_indices)
        return false
    end

    # Build adjacency list and find active vertices
    adj_list = Dict(i => [] for i in 1:8)
    active_vertices = Set{Int}()

    for edge_idx in active_edge_indices
        u, v = edges_array[edge_idx, :]
        push!(adj_list[u], v)
        push!(adj_list[v], u)
        push!(active_vertices, u)
        push!(active_vertices, v)
    end

    # Perform a Depth-First Search (DFS)
    start_node = first(active_vertices)
    stack = [start_node]
    visited = Set{Int}()

    while !isempty(stack)
        node = pop!(stack)
        if !(node in visited)
            push!(visited, node)
            for neighbor in adj_list[node]
                push!(stack, neighbor)
            end
        end
    end

    return length(visited) == length(active_vertices)
end

# --- VISUALIZATION FUNCTION ---
function plot_all_incomplete_cubes(cubes_to_plot; active_sides=nothing, require_connected=false)
    # 1. Filter cubes based on criteria
    title_parts = ["Rotationally Unique Cubes"]
    filtered_cubes = cubes_to_plot

    if active_sides !== nothing
        push!(title_parts, "with $active_sides Edges")
        filter!(c -> sum(c) == active_sides, filtered_cubes)
    end

    if require_connected
        pushfirst!(title_parts, "Connected")
        filter!(c -> is_connected(c, edges), filtered_cubes)
    end

    title_str = join(title_parts, " ")

    if isempty(filtered_cubes)
        println("No cubes found matching the criteria.")
        return
    end

    # 2. Generate coordinates for each cube with an offset
    x_coords, y_coords, z_coords = [], [], []
    x_offset = 0.0
    spacing = 2.0

    for visibility_vector in filtered_cubes
        for (i, is_visible) in enumerate(visibility_vector)
            if is_visible == 1
                v_indices = edges[i, :]
                start_v = vertices[v_indices[1], :]
                end_v = vertices[v_indices[2], :]

                append!(x_coords, [start_v[1] + x_offset, end_v[1] + x_offset, nothing])
                append!(y_coords, [start_v[2], end_v[2], nothing])
                append!(z_coords, [start_v[3], end_v[3], nothing])
            end
        end
        x_offset += spacing
    end

    # 3. Create and configure the plot
    trace = scatter3d(
        x=x_coords, y=y_coords, z=z_coords,
        mode="lines",
        line=attr(color="#22d3ee", width=5),
        hoverinfo="none"
    )

    x_range = x_offset - spacing + 2.0
    y_range = 2.0
    z_range = 2.0

    layout = Layout(
        title=title_str,
        paper_bgcolor="#111827",
        plot_bgcolor="#111827",
        margin=attr(l=0, r=0, b=0, t=40),
        scene=attr(
            xaxis=attr(visible=false, range=[-0.5, x_offset - spacing + 1.5]),
            yaxis=attr(visible=false, range=[-0.5, 1.5]),
            zaxis=attr(visible=false, range=[-0.5, 1.5]),
            aspectmode="manual",
            aspectratio=attr(x=x_range, y=y_range, z=z_range)
        )
    )

    display(plot(trace, layout))
end

# --- MAIN EXECUTION ---
function main()
    cube_data_path = "new_cubes.npy"

    if !isfile(cube_data_path)
        println("Error: Data file not found at '$cube_data_path'")
        println("Please run 'compute_cubes.jl' first to generate this file.")
        return
    end

    all_unique_cubes_matrix = npzread(cube_data_path)
    # Convert matrix to a vector of vectors for easier iteration
    all_unique_cubes = [all_unique_cubes_matrix[i, :] for i in 1:size(all_unique_cubes_matrix, 1)]

    println("Loaded $(length(all_unique_cubes)) unique cube configurations.")

    active_sides_to_plot = 7

    println("\nPlotting all CONNECTED unique cubes with $active_sides_to_plot active edges...")
    plot_all_incomplete_cubes(
        all_unique_cubes,
        active_sides=active_sides_to_plot,
        require_connected=true
    )
end

# Run the main function
main()

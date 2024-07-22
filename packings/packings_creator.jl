# The purpose of this script is to create data files
# for LAMMPS granular packing simulations. Using Julia,
# you can call this function using 
# MakeLammpsICs(N, W) in the command line where you
# determine N and W as numbers.
#
# W = width of square particle top and bottom boundaries
# defined by the number of particles EXCLUDING the particles
# at the boundary. For example, a W = 5 will create a 4x4 wall.
#
# N = number of particles that are between the walls.
#
# Example, a "MakeLammpsICs(1000, 5)" call will create a 
# packing with 4x4x2=32 particles on the top and bottom walls
# with 1000 particles between them for a total of 1032 
# particles.
####################################################
using Random
using Plots
plotlyjs()

function create_packings_3D(N,W, diameter_average, diameter_spread, bi_disperse)
    

    if bi_disperse
        diameter_large = diameter_average + diameter_spread
        diameter_small = diameter_average
    end

    # Calulate the heigh to of the simulation box and scale to the size of the largest particle
    x_limit_upper = W*diameter_large
    y_limit_upper = W*diameter_large
    z_limit_upper = N/(W^2)*diameter_large

    ##  Create a rank 3 tensor that is x_limit_upper by y_limit_upper by z_limit_upper where each point is the x,y,z coordinate that is evenly spaced by diameter_large 
    # Generate ranges for x, y, and z coordinates so they don't extend past the boundary
    x_range = collect(range(diameter_large/2, x_limit_upper-diameter_large/2, step=diameter_large)) 
    y_range = collect(range(diameter_large/2, y_limit_upper-diameter_large/2, step=diameter_large))
    z_range = collect(range(diameter_large/2, z_limit_upper-diameter_large/2, step=diameter_large))

    # Create 3D grid of coordinates
    x_coords = [x for x in x_range, y in y_range, z in z_range] # creates length(z_range) matricies that have dimensions length(x) by length(y) where each element is the x-position
    y_coords = [y for x in x_range, y in y_range, z in z_range]
    z_coords = [z for x in x_range, y in y_range, z in z_range]

    # concatanate into a matrix where each row is the x,y,z coordinate of each particle has 
    coordinate_list = hcat(reshape(x_coords, :, 1), reshape(y_coords, :, 1), reshape(z_coords, :, 1))

    
    # Readjust the number of particles to those that actually fit in the box
    N_actual = size(coordinate_list,1)


    ## Take the new N and then randomly assign a partcile diamter to reach between diameter_large and diameter_small
    particle_diameters = diameter_small .+ (diameter_large - diameter_small) .* rand(N_actual)


    # Create the 3D plot
    x = coordinate_list[:, 1]
    y = coordinate_list[:, 2]
    z = coordinate_list[:, 3]
    marker_sizes = particle_diameters * 10  # Adjust the scaling factor if needed for better visualization

    # Generate hover texts including particle number and coordinates
    hover_texts = ["Particle $i: (x, y, z) = ($(x[i]), $(y[i]), $(z[i]))" for i in 1:N_actual]
    scatter3d(x, y, z, markersize=marker_sizes, label="", xlabel="X", ylabel="Y", zlabel="Z", title="3D Scatter Plot of Particles", hover=hover_texts)

    min_z = minimum(coordinate_list[:,3])
    max_z = maximum(coordinate_list[:,3])
    indicies_bottom_layer = findall(z -> z == min_z, coordinate_list[:, 3])
    indicies_top_layer = findall(z -> z == max_z, coordinate_list[:, 3])

    # Write data to input_file.lammps
    open("input_file.lammps", "w") do file

    # Write header information
    write(file, "LAMMPS 3d granular data file written by packigns_creator.jl\n\n")
    write(file, "$N_actual atoms\n0 bonds\n0 angles\n0 dihedrals\n0 impropers\n\n")
    write(file, "2 atom types\n\n")
    write(file, "0.0 $x_limit_upper xlo xhi\n")
    write(file, "0.0 $y_limit_upper ylo yhi\n")
    write(file, "0.0 $z_limit_upper zlo zhi\n\n")
    
    # Write atoms section
    write(file, "Atoms\n\n")
        for i in indicies_bottom_layer[1]:indicies_bottom_layer[end]
            # Format is atom_ID, atom_type, diameter, rho, x, y, z
            write(file, "$i 2 $(particle_diameters[i]) 1 $(x[i]) $(y[i]) $(z[i]) \n")
        end

        for i in indicies_bottom_layer[end]+1:indicies_top_layer[1]-1
            # Format is atom_ID, atom_type, diameter, rho, x, y, z
            write(file, "$i 1 $(particle_diameters[i]) 1 $(x[i]) $(y[i]) $(z[i]) \n")
        end

        for i in indicies_top_layer[1]:indicies_top_layer[end]
            # Format is atom_ID, atom_type, diameter, rho, x, y, z
            write(file, "$i 3 $(particle_diameters[i]) 1 $(x[i]) $(y[i]) $(z[i]) \n")
        end
    end
end



N = 100
W = 5
diameter_average = 1
diameter_spread = .4
bi_disperse = true


N= create_packings_3D(N,W, diameter_average, diameter_spread, bi_disperse)
println("$N")

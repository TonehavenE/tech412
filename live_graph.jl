using LinearAlgebra, GLMakie
using GLMakie.GLFW
using GLMakie: to_native

"""
    rotation(θ)

return 2 x 2 rotation matrix that rotates plane by angle θ
"""
function rotation(θ)
    c, s = cos(θ), sin(θ)
    return [c -s; s c]
end

"""
    dihedralgroup(n, flip=true)

Generate symmetries of the n-gon that rotate and flip the n-gon in the plane.
The return value is an array of n or 2n matrices representing the elements of the group.
"""
function dihedralgroup(n, flip=true)    
    S = [-1 0; 0 1]                 # a reflection about y axis
    I = [1.0 0.0; 0.0 1.0]          # the identity
    Dn = fill(I, flip ? 2n : n) # allocate an array of 2n or n matrices
    
    for k=1:n
        @views Dn[k] = rotation(2(k-1)π/n) # set Dn[k] to rotation by θ = 2(k-1)π/n
        if flip
            Dn[k+n] = S*Dn[k]       # set Dn[k+n] to reflection of Dn[k]
        end
    end
    Dn
end

"""
    symmetrize(X, G)

Symmetrize a set of data points X by symmetry group G. The return value
is a matrix containing all columns of X mapped by all matrices in G
"""
function symmetrize(X, G)
    m, nX = size(X)
    nG = length(G)
    GX = similar(X, m, nX * nG)  # Preallocate memory

    for j in 1:nX
        @views Xj = X[:, j]  # Avoid unnecessary indexing
        for k in 1:nG
            GX[:, (j-1) * nG + k] = G[k] * Xj
        end
    end
    return GX
end


"""
    f(x, X, a=1, k=1)

return 1/N sum_j cos(k|x-xj|) exp(-a|x-xj|^2) where 
  x is a 2d vector 
  xj is the jth column of 2 x N matrix X
"""
function f(x, X, a=1, k=1)
	norms = vec(norm.(eachcol(X .- x)))
	return sum(cos.(k .* norms) .* exp.(-a .* norms.^2))/ size(X, 2)
end

using Base.Threads

function generate_zgrid(Xsymm, a, k, x1grid, x2grid)
	N1 = length(x1grid)
	N2 = length(x2grid)
	zgrid = zeros(N1, N2)

	@threads for i in 1:N1
		for j in 1:N2
			x = [x1grid[i]; x2grid[j]]
			zgrid[i, j] = f(x, Xsymm, a, k)
		end
	end
	return zgrid
end

function compute_pattern(X, a, k, w, n, flip)	
    G = dihedralgroup(n, flip)
    Xsymm = symmetrize(X, G)

    x1grid = range(-w, w, length=100)
    x2grid = range(-w, w, length=100)

    # Generate zgrid based on the input data and current parameters
    zgrid = generate_zgrid(Xsymm, a, k, x1grid, x2grid)
    zscale = maximum(abs.(zgrid))

    return zgrid / zscale
end

# Function to create animation with sliders
function observable_animation()
    w = 400
    fps = 60
    polygon = 5
    a0 = 3
    k0 = 5
    s = 3
    Npts = 3
    width = 3
    flip = true
    levels = -1:.3:1

    X = randn(2, Npts)
    Xt = copy(X)
    Nx = size(X, 2)
    # time = Observable(0.0)  # Observable for time
    dt = pi / 512

    # Movement functions as observables
    a(a0, t) = a0*(1 .- 2/3*cos.(t))
    k(k0, t) = k0*(1 .- 2/3*cos.(2*t))

    # Angular velocities for rotation
    ω = (rand(Int, Npts) .% 4 .+ 1).* (-1).^(rand(Int, Npts) .% 2)
    if !flip
        ω .= ω .- sum(ω)/Npts  # remove mean rotation
    end

    R = rotation.(ω * dt)  # Precompute rotation matrices

    # Observable for the zgrid data
    zgrid = Observable(fill(0.0, 100, 100))  # Placeholder, will be updated with actual data

    # Create figure and axis
    fig = Figure(size=(1920, 1080))
    ax = Axis(fig[1, 1], aspect=1, xgridvisible=false, ygridvisible=false)

    # Enable quitting with "q" key
    glfw_window = to_native(display(fig))
    on(events(fig).keyboardbutton) do event
        if event.key == Keyboard.q
            GLFW.SetWindowShouldClose(glfw_window, true)
        end
    end

    # Sliders for interactive parameters
    slider_grid = SliderGrid(
        fig[2, 1], 
        (label="a0", range=0:0.01:10, startvalue=3),
        (label="k0", range=0:0.01:10, startvalue=5),
        (label="Polygon", range=0:1:10, startvalue=5),
    )
    sliderobservables = [s.value for s in slider_grid.sliders]

    # Colormap menu
    menu = Menu(fig[3, 1], options=["viridis", "plasma", "inferno", "magma", "coolwarm", "turbo"], default="viridis")
    selected_colormap = menu.selection
    contourf!(ax, range(-w, w, length=100), range(-w, w, length=100), zgrid, colormap=selected_colormap[], levels=levels)
    display(fig)

    # Update the observable zgrid whenever the time or parameters change
	nframes = 100
	time = 0
	for i = 1:nframes
		if GLFW.WindowShouldClose(glfw_window)
			break
		end
        a0, k0, polygon = (s[] for s in sliderobservables)
        # Calculate new pattern
        new_zgrid = compute_pattern(Xt, a(a0, time), k(k0, time), width, polygon, flip)

        # Update zgrid observable
        zgrid[] = new_zgrid

		sleep(1/fps)
		time += 1/fps
	end
end

# animation_with_sliders()
observable_animation()

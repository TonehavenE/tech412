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
    s = 0.0
    N = size(X, 2)
    for j in 1:N
        r = norm(x-X[:,j])
        s += cos(k*r)*exp(-a*(r^2))
    end
    s/N
end
# function f(x, X, a=1, k=1)
#     r = @views norm.(eachcol(X) .- x)  # Compute all norms at once
#     return sum(cos.(k .* r) .* exp.(-a .* r.^2)) / size(X, 2)
# end

"""
    plotpattern(n, flip, X, a, k, width, levels, colormap)    
     
Generate a symmetric pattern based on 
  n : use symmetry group of n-gon
  flip : boolean, use / don't use reflection symmetries
  X : 2 x N matrix of data points, each column of X is a point in the plane 
  a : scale of blobs, exp(-a r^2)
  k : scale of oscillations, cos(k r^2)
  width : width and height of plot axes
  levels : number or values of contour levels
  colormap : colormap for contour plot     
"""
function plotpattern!(ax, n, flip, X, a, k, width, levels, colormap)
    w = width
    # define our groups and points
    G = dihedralgroup(n, flip)
    Xsymm = symmetrize(X, G)

    # evaluate f(x, Xsymm, a, k) over a grid of points x=[x1;x2]
    
    x1grid = range(-w, w, length=100)
    x2grid = range(-w, w, length=100)
    # symmetrize the data points X with the symmetries of the n-gon
    zgrid = [f([x1;x2], Xsymm, a, k) for x2 in x2grid, x1 in x1grid]
    zscale = maximum(abs.(zgrid))
    # make a contour plot of zgrid = f(x, Xsymm, a, k)
    empty!(ax)
    contourf!(ax, x1grid, x2grid, zgrid/zscale, colormap=colormap, levels=levels)
end

"""
Create an animation with sliders.
"""
function animation_with_sliders()
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
    time = 0
    dt = pi / 512

	# movement functions
	a(a0, t) = a0*(1 .- 2/3*cos.(t))
	k(k0, t) = k0*(1 .- 2/3*cos.(2*t));

	ω = (rand(Int,Npts) .% 4 .+ 1).* (-1).^(rand(Int,Npts) .% 2)
	if !flip
	    ω = ω .- sum(ω)/Npts      # remove mean rotation
	end

    R = rotation.(ω * dt)  # Precompute rotation matrices

    # Create figure
    fig = Figure(size=(1920, 1080))
    ax = Axis(fig[1, 1], aspect=1, xgridvisible=false, ygridvisible=false)
    display(fig)

    # Enable quitting
    glfw_window = to_native(display(fig))
    on(events(fig).keyboardbutton) do event
        if event.key == Keyboard.q
            GLFW.SetWindowShouldClose(glfw_window, true)
        end
    end

    # Sliders
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


    while !GLFW.WindowShouldClose(glfw_window)
		# println("plotting new frame")
        a0, k0, polygon = (s[] for s in sliderobservables)
			plotpattern!(ax, polygon, flip, Xt, a(a0, time), k(k0, time), width, levels, selected_colormap[])
        time += 1 / fps
        sleep(1 / fps)

        # @views Xt .= R .* Xt  # In-place rotation
		# Xt .= R .* Xt
		for j=1:Nx
			Xt[:, j] = R[j]*Xt[:, j]
		end
    end

    return fig
end


animation_with_sliders()

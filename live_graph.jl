using LinearAlgebra, GLMakie
using GLMakie.GLFW
using GLMakie: to_native

"""
    rotation(θ)

return 2 x 2 rotation matrix that rotates plane by angle θ
"""
rotation(θ) = [cos(θ) -sin(θ); sin(θ) cos(θ)]

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
        Dn[k] = rotation(2(k-1)π/n) # set Dn[k] to rotation by θ = 2(k-1)π/n
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
    m,nX = size(X)  # nX is number of data points
    nG = length(G)  # nG is number elements in group
    
    GX = fill(0.0, m, nX*nG) # allocate a matrix for G applied to X
    
    for j in 1:nX      # for each datapoint in X...
        for k in 1:nG  # ...and for each matrix in the group...
            GX[:, (j-1)*nG + k] = G[k]*X[:,j] # ...map the jth datapoint by the kth matrix
        end
    end
    GX
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
    nframes = 240
    time = 0
    polygon = 5
    a0 = 3    # scale of blobs (larger a, narrower blobs)
    k0 = 5    # scale of ripples (larger k, more rapid ripples)
    s = 3     # scale of data points (larger s, further spread out)
	Npts = 3
	width = 3
	flip = true
	levels = -1:.3:1  # number or values of contour levels

	X = randn(2, Npts)
	dt = pi/512
	t = 0:dt:2pi

	# some choices for rotation rates
	#ω = 1.5*(2*rand(Npts).-1) # random rotation rates for data points 
	#ω = ω .- sum(ω)/Npts      # remove mean rotation
	#ω = rand(Int, Npts) .% 6 .+ 1
	ω = (rand(Int,Npts) .% 4 .+ 1).* (-1).^(rand(Int,Npts) .% 2)
	if !flip
	    ω = ω .- sum(ω)/Npts      # remove mean rotation
	end

	a(a0, t) = a0*(1 .- 2/3*cos.(t))
	k(k0, t) = k0*(1 .- 2/3*cos.(2*t));

	R = rotation.(ω*dt)   # vector of incremental rotation matrices, one for each data point
 
    Nx = size(X,2)
    Xt = copy(X)
    
	# create figure
	fig = Figure(size=(1920, 1080))
    ax = Axis(fig[1, 1], aspect = 1, xgridvisible = false, ygridvisible = false)
    display(fig)
	
	# enable quitting
	glfw_window = to_native(display(fig))

	on(events(fig).keyboardbutton) do event
    	if event.key == Keyboard.q
      		GLFW.SetWindowShouldClose(glfw_window, true) # this will close the window after all callbacks are finished
		end
    end

    # Define slider grid
    slider_grid = SliderGrid(
        fig[2, 1], 
        (label = "a0", range = 0:0.01:10, startvalue = 3),
        (label = "k0", range = 0:0.01:10, startvalue = 5),
        (label = "Polygon", range = 0:1:10, startvalue = 5),
    )
    sliderobservables = [s.value for s in slider_grid.sliders]

    # Define colormap menu
    colormap_options = ["viridis", "plasma", "inferno", "magma", "coolwarm", "turbo"]
    menu = Menu(fig[3, 1], options=colormap_options, default="viridis")
    selected_colormap = menu.selection  # Observable storing the selected colormap


    # for i = 1:nframes
	while true
        a0, k0, polygon = sliderobservables[1][], sliderobservables[2][], sliderobservables[3][]
        plotpattern!(ax, polygon, flip, Xt, a(a0, time), k(k0, time), width, levels, selected_colormap[])
        time += 1/fps
        sleep(1/fps)

		for j=1:Nx
			Xt[:, j] = R[j]*Xt[:, j]
		end
    end

    return fig
end

animation_with_sliders()

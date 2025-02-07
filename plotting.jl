using LinearAlgebra, GLMakie

println("Initalized")
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


function plotpattern(n, flip, X, a, k, width, levels, colormap) 
    # symmetrize the data points X with the symmetries of the n-gon
    G = dihedralgroup(n, flip)
    Xsymm = symmetrize(X, G)

    # evaluate f(x, Xsymm, a, k) over a grid of points x=[x1;x2]
    w = width 
    x1grid = range(-w, w, length=100)
    x2grid = range(-w, w, length=100)
    zgrid = [f([x1;x2], Xsymm, a, k) for x2 in x2grid, x1 in x1grid]

    # make a contour plot of zgrid = f(x, Xsymm, a, k)
    zscale = maximum(abs.(zgrid))
    fig = Figure(size=(400, 400))
    ax = Axis(fig[1, 1], aspect = 1)
    contourf!(ax, x1grid, x2grid, zgrid/zscale, colormap=colormap, levels=levels)
    xlims!(ax, -w, w)
    ylims!(ax, -w, w)
	
	return fig
end

# pattern parameters
ngon = 5  # n-gon dihedral group
Npts = 3  # number of random points
a0 = 3    # scale of blobs (larger a, narrower blobs)
k0 = 5    # scale of ripples (larger k, more rapid ripples)
s = 3     # scale of data points (larger s, further spread out)
width = 3 # width of plot: -width < x < width, -width < y < width

flip = true  # do or don't include mirror symmetry 
speed = 16   # if animation is too fast, reduce this value
levels = -1:.3:1  # number or values of contour levels
colormap =  :ice # :Set1 # :seaborn_colorblind6 # :Set1_4 #:PuRd_5  #:seaborn_colorblind6  #:Paired_8 (search on "Julia Plots colormaps" to find other color palettes)

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

# time variation of blob scale a(t) and ripple scale k(t)
a(t) = a0*(1 .- 2/3*cos.(t))
k(t) = k0*(1 .- 2/3*cos.(2*t));

fig = plotpattern(ngon, flip, X, a(t[1]), k(t[1]), width, levels, colormap)
display(fig)


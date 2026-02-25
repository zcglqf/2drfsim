# make_sphere2d_phantom1024.jl
# How to use: in Julia REPL, include("make_sphere2d_phantom1024.jl"), include full path if not in current directory.
# and call: phantom = make_phantom_2d(; Nx=256, Ny=256, FOVx=0.24, FOVy=0.24)

using KomaMRI

"""
Create a 2D KomaMRI Phantom on a centered (x,y) grid.

Units
  FOVx, FOVy in meters
  x,y,z in meters
  T1,T2 in seconds

Returns
  phantom::Phantom
"""
function make_phantom_2d(; Nx::Int=256, Ny::Int=256,
                          FOVx::Float64=0.256, FOVy::Float64=0.256,
                          rho_val::Float64=1.0)

    dx = FOVx / Nx
    dy = FOVy / Ny

    # Centered grid, sample at voxel centers
    x_lin = range(-FOVx/2 + dx/2, FOVx/2 - dx/2, length=Nx)
    y_lin = range(-FOVy/2 + dy/2, FOVy/2 - dy/2, length=Ny)

    # Flattened coordinate vectors (length = Nx*Ny)
    X = repeat(collect(x_lin), inner=Ny)
    Y = repeat(collect(y_lin), outer=Nx)
    Z = zeros(length(X))

    Nspins = length(X)
    ρ  = fill(rho_val, Nspins)

    phantom = Phantom(name="HR2D",x=X, y=Y, z=Z, ρ=ρ)
    return phantom
end

"""
Optional helper: make a circular object with different T1/T2/rho inside radius r (meters).
Call this after make_phantom_2d.
"""
function apply_circle!(phantom; r::Float64=0.05,
                       T1_in::Float64=1400e-3, T2_in::Float64=100e-3, rho_in::Float64=1.0,
                       T1_out::Float64=0.0,    T2_out::Float64=0.0,   rho_out::Float64=0.0)

    x = getproperty(phantom, :x)
    y = getproperty(phantom, :y)

    rr = x .* x .+ y .* y
    inside = rr .<= r^2

    phantom.T1 .= T1_out
    phantom.T2 .= T2_out
    getproperty(phantom, Symbol("ρ")) .= rho_out

    phantom.T1[inside] .= T1_in
    phantom.T2[inside] .= T2_in
    getproperty(phantom, Symbol("ρ"))[inside] .= rho_in

    return phantom
end
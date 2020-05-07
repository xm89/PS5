include("PhasePortraitV2.jl")

# Function for a dual repression system without cooperativity
# x1: range of x1 values (i.e. A values)
# x2: range of x2 values (i.e. R values)
# We use `@.` to apply the calculations across all rows.
# Note that model parameters are specified within the function
# Returns computed (dx1/dt, dx2/dt) over the range of (x1, x2)
function toggleMono(x1, x2)
    alpha = 10.0               #degradation rate const. for repressor 1
    n = 2.0               #degradation rate const. for repressor 2

    u = @. -x1 + alpha/(1+x2^n)
    v = @. -x2 + alpha/(1+x1^n)   #eqn 12

    return (u,v)
end

#Range of x1, x2 values
x1range = (0,15,15)          #Has the form (min, max, num points)
x2range = (0,15,15)          #Has the form (min, max, num points)
x₀ = ([1.0, 10.0],)  #initial state vectors; a common must be included after the first array
tspan=(0.0,30.0)             #time span

#Call the phaseplot functon to construct the phase portrait
phaseplot(toggleMono, x1range, x2range, xinit=(), t=tspan, clines=true,
        norm=true, scale=0.5)

using DrWatson
@quickactivate "project"
include(srcdir("SIRPetri.jl"))
using .SIRPetri
using Plots

β = 0.3
γ = 0.1
tmax = 100.0

net, u0, _ = build_sir_network(β, γ)
df = simulate_deterministic(net, u0, (0.0, tmax), saveat = 0.5, rates = [β, γ])

anim = Animation()
for i in 1:5:length(df.time)
    p = bar(
        ["S", "I", "R"],
        [df.S[i], df.I[i], df.R[i]],
        ylims = (0, 1000),
        title = "Time = $(round(df.time[i], digits=1))",
        ylabel = "Population",
        label = "",
        color = [:green, :red, :blue]
    )
    frame(anim, p)
end

gif(anim, plotsdir("sir_animation.gif"), fps = 5)

println("Анимация сохранена в plots/sir_animation.gif")

using DrWatson
@quickactivate "project"
using DataFrames, CSV, Plots

df_det = CSV.read(datadir("sir_det.csv"), DataFrame)
df_stoch = CSV.read(datadir("sir_stoch.csv"), DataFrame)
df_scan = CSV.read(datadir("sir_scan.csv"), DataFrame)

n_points = min(length(df_det.time), length(df_stoch.time))
p1 = plot(
    df_det.time[1:n_points],
    [df_det.I[1:n_points] df_stoch.I[1:n_points]],
    label = ["Deterministic I" "Stochastic I"],
    xlabel = "Time",
    ylabel = "Infected",
    title = "Comparison of Deterministic and Stochastic Dynamics",
    linewidth = 2,
)
savefig(plotsdir("comparison.png"))

p2 = plot(
    df_scan.β,
    df_scan.peak_I,
    marker = :circle,
    linewidth = 2,
    xlabel = "β (infection rate)",
    ylabel = "Peak Infected",
    title = "Sensitivity Analysis: Peak I vs β",
)
savefig(plotsdir("sensitivity.png"))

println("Отчётные графики сохранены в plots/")

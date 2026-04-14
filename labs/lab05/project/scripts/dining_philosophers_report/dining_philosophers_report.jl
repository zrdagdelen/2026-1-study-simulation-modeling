using DrWatson
@quickactivate "project"

using DataFrames, CSV, Plots

df_classic = CSV.read(datadir("dining_classic.csv"), DataFrame)
df_arbiter = CSV.read(datadir("dining_arbiter.csv"), DataFrame)

N = 5

p1 = plot(
    xlabel = "Время",
    ylabel = "Ест (1/0)",
    title = "Классическая сеть (без арбитра)",
    legend = :topright,
    legendtitle = "Философ",
)

for i = 1:N
    col = Symbol("Eat_$i")
    plot!(p1, df_classic.time, df_classic[:, col], label = "Ф $i")
end

p2 = plot(
    xlabel = "Время",
    ylabel = "Ест (1/0)",
    title = "Сеть с арбитром",
    legend = :topright,
    legendtitle = "Философ",
)

for i = 1:N
    col = Symbol("Eat_$i")
    plot!(p2, df_arbiter.time, df_arbiter[:, col], label = "Ф $i")
end

p_final = plot(p1, p2, layout = (2, 1), size = (800, 600))
savefig(plotsdir("final_report.png"))

println("Отчёт сохранён в plots/final_report.png")

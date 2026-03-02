using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using DifferentialEquations
using DataFrames
using StatsPlots
using Plots
using LaTeXStrings
using Statistics

function lotka_volterra!(du, u, p, t)
    x, y = u
    α, β, δ, γ = p

    du[1] = α*x - β*x*y
    du[2] = δ*x*y - γ*y
    nothing
end

u0 = [40.0, 9.0]
tspan = (0.0, 200.0)

p_base = [0.1, 0.02, 0.01, 0.3]

function run_lv(p; tspan_input=tspan)
    prob = ODEProblem(lotka_volterra!, u0, tspan_input, p)
    sol = solve(prob, dt=0.1, saveat=0.5)

    df = DataFrame()
    df.t = sol.t
    df.prey = [u[1] for u in sol.u]
    df.predator = [u[2] for u in sol.u]

    x_star = p[4] / p[3]
    y_star = p[1] / p[2]

    return df, x_star, y_star
end

df_base, x_star_base, y_star_base = run_lv(p_base)

println("="^60)
println("БАЗОВЫЙ СЦЕНАРИЙ")
println("="^60)
println("Параметры: α=$(p_base[1]), β=$(p_base[2]), δ=$(p_base[3]), γ=$(p_base[4])")
println("Стационарная точка: x* = $(round(x_star_base, digits=2)) (жертвы), y* = $(round(y_star_base, digits=2)) (хищники)")
println("Средняя численность жертв: $(round(mean(df_base.prey), digits=2))")
println("Средняя численность хищников: $(round(mean(df_base.predator), digits=2))")

println("\n" * "="^60)
println("ИССЛЕДОВАНИЕ 1: Влияние скорости размножения жертв (α)")
println("="^60)

α_values = [0.05, 0.1, 0.2, 0.3]
β_fixed, δ_fixed, γ_fixed = p_base[2], p_base[3], p_base[4]

plt_α = plot(layout=(2,1), size=(800,600), legend=:topright)

for α in α_values
    p = [α, β_fixed, δ_fixed, γ_fixed]
    df, x_star, y_star = run_lv(p)

    plot!(plt_α[1], df.t, df.prey, label="α=$α", linewidth=2)
    plot!(plt_α[2], df.t, df.predator, label="α=$α", linewidth=2)

    println("α = $α: стац. точка: x*=$(round(x_star, digits=1)), y*=$(round(y_star, digits=2))")
end

plot!(plt_α[1], title="Динамика популяции жертв", ylabel="Численность жертв", grid=true)
plot!(plt_α[2], title="Динамика популяции хищников", xlabel="Время", ylabel="Численность хищников", grid=true)

println("\n" * "="^60)
println("ИССЛЕДОВАНИЕ 2: Влияние эффективности охоты (β)")
println("="^60)

β_values = [0.01, 0.02, 0.03, 0.04]
α_fixed, δ_fixed, γ_fixed = p_base[1], p_base[3], p_base[4]

plt_β = plot(layout=(2,1), size=(800,600), legend=:topright)

for β in β_values
    p = [α_fixed, β, δ_fixed, γ_fixed]
    df, x_star, y_star = run_lv(p)

    plot!(plt_β[1], df.t, df.prey, label="β=$β", linewidth=2)
    plot!(plt_β[2], df.t, df.predator, label="β=$β", linewidth=2)

    println("β = $β: стац. точка: x*=$(round(x_star, digits=1)), y*=$(round(y_star, digits=2))")
end

plot!(plt_β[1], title="Динамика жертв", ylabel="Численность жертв", grid=true)
plot!(plt_β[2], title="Динамика хищников", xlabel="Время", ylabel="Численность хищников", grid=true)

println("\n" * "="^60)
println("ИССЛЕДОВАНИЕ 3: Влияние конверсии биомассы (δ)")
println("="^60)

δ_values = [0.005, 0.01, 0.02, 0.03]
α_fixed, β_fixed, γ_fixed = p_base[1], p_base[2], p_base[4]

plt_δ = plot(layout=(2,1), size=(800,600), legend=:topright)

for δ in δ_values
    p = [α_fixed, β_fixed, δ, γ_fixed]
    df, x_star, y_star = run_lv(p)

    plot!(plt_δ[1], df.t, df.prey, label="δ=$δ", linewidth=2)
    plot!(plt_δ[2], df.t, df.predator, label="δ=$δ", linewidth=2)

    println("δ = $δ: стац. точка: x*=$(round(x_star, digits=1)), y*=$(round(y_star, digits=2))")
end

plot!(plt_δ[1], title="Динамика жертв", ylabel="Численность жертв", grid=true)
plot!(plt_δ[2], title="Динамика хищников", xlabel="Время", ylabel="Численность хищников", grid=true)

println("\n" * "="^60)
println("ИССЛЕДОВАНИЕ 4: Влияние смертности хищников (γ)")
println("="^60)

γ_values = [0.1, 0.3, 0.5, 0.7]
α_fixed, β_fixed, δ_fixed = p_base[1], p_base[2], p_base[3]

plt_γ = plot(layout=(2,1), size=(800,600), legend=:topright)

for γ in γ_values
    p = [α_fixed, β_fixed, δ_fixed, γ]
    df, x_star, y_star = run_lv(p)

    plot!(plt_γ[1], df.t, df.prey, label="γ=$γ", linewidth=2)
    plot!(plt_γ[2], df.t, df.predator, label="γ=$γ", linewidth=2)

    println("γ = $γ: стац. точка: x*=$(round(x_star, digits=1)), y*=$(round(y_star, digits=2))")
end

plot!(plt_γ[1], title="Динамика жертв", ylabel="Численность жертв", grid=true)
plot!(plt_γ[2], title="Динамика хищников", xlabel="Время", ylabel="Численность хищников", grid=true)

println("\n" * "="^60)
println("ФАЗОВЫЕ ПОРТРЕТЫ")
println("="^60)

plt_phase = plot(layout=(2,2), size=(900,800))

params_phase = [
    ([0.1, 0.02, 0.01, 0.3], "Базовый"),
    ([0.2, 0.02, 0.01, 0.3], "Высокий α"),
    ([0.1, 0.04, 0.01, 0.3], "Высокий β"),
    ([0.1, 0.02, 0.02, 0.3], "Высокий δ")
]

for (i, (p, title)) in enumerate(params_phase)
    df, x_star, y_star = run_lv(p; tspan_input=(0.0, 300.0))

    plot!(plt_phase[i], df.prey, df.predator, linewidth=2, label="Траектория", color=:blue)
    scatter!(plt_phase[i], [x_star], [y_star], color=:red, markersize=8, label="Стац. точка")
    plot!(plt_phase[i], title=title, xlabel="Жертвы", ylabel="Хищники", grid=true, legend=:topright)

    println("$title: стац. точка (x=$(round(x_star, digits=1)), y=$(round(y_star, digits=2)))")
end

plots_path = joinpath(@__DIR__, "..", "plots", "lv_literate")
mkpath(plots_path)

savefig(plt_α, joinpath(plots_path, "lv_alpha_study.png"))
savefig(plt_β, joinpath(plots_path, "lv_beta_study.png"))
savefig(plt_δ, joinpath(plots_path, "lv_delta_study.png"))
savefig(plt_γ, joinpath(plots_path, "lv_gamma_study.png"))
savefig(plt_phase, joinpath(plots_path, "lv_phase_portraits.png"))

println("\nВсе графики сохранены в: ", plots_path)

println("\n" * "="^60)
println("СВОДНАЯ ТАБЛИЦА РЕЗУЛЬТАТОВ")
println("="^60)

println("\nВлияние параметров на стационарную точку:")
println("┌──────────┬───────────────┬───────────────┐")
println("│ Параметр │ Влияние на x* │ Влияние на y* │")
println("├──────────┼───────────────┼───────────────┤")
println("│    α ↑   │      —        │      ↑        │")
println("│    β ↑   │      —        │      ↓        │")
println("│    δ ↑   │      ↓        │      —        │")
println("│    γ ↑   │      ↑        │      —        │")
println("└──────────┴───────────────┴───────────────┘")

println("\n" * "="^60)
println("Исследование модели Лотки-Вольтерры завершено!")
println("="^60)

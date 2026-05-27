#!/usr/bin/env julia

using DrWatson
@quickactivate "project"

include(srcdir("sir_model.jl"))

using Random      # Генерация случайных чисел и фиксация seed
using StatsPlots  # Построение графиков
using BenchmarkTools # Оценка производительности
using CSV         # Сохранение результатов в CSV
using Dates       # Работа с датами для именования файлов

tmax = 40.0

uθ = [990, 10, 0]  # S, I, R

p_default = [0.05, 10.0, 0.25]  # β, c, γ

Random.seed!(1234)

println("="^60)
println("ДИСКРЕТНО-СОБЫТИЙНАЯ SIR МОДЕЛЬ")
println("="^60)

println("\n1. Базовый запуск модели...")

des_model = MakeSIRModel(uθ, p_default)

activate(des_model)

sir_run(des_model, tmax)

data_des = out(des_model)

plot(data_des.t, [data_des.S data_des.I data_des.R],
     labels=["S" "I" "R"],
     xlabel="Время",
     ylabel="Численность",
     title="Дискретно-событийная SIR модель (β=0.05, c=10.0, γ=0.25)",
     linewidth=2)

savefig(plotsdir("sir_des.png"))

println("   График сохранён: plots/sir_des.png")

println("\n2. Анализ чувствительности к параметрам...")

betas = [0.03, 0.05, 0.07]
results_beta = []

println("\n   Варьирование β (c=10.0, γ=0.25):")
for β in betas

    Random.seed!(1234)
    p = [β, 10.0, 0.25]
    m = MakeSIRModel(uθ, p)
    activate(m)
    sir_run(m, tmax)
    data = out(m)

    peak_I = maximum(data.I)                        # Пик инфицированных
    peak_time = data.t[argmax(data.I)]              # Время пика
    final_R = data.R[end]                           # Итоговое число переболевших

    push!(results_beta, (β=β, peak_I=peak_I, peak_time=peak_time, final_R=final_R))
    println("   β=$β: пик I=$(round(peak_I)), время=$(round(peak_time, digits=2)), переболело=$(round(final_R))")
end

cs = [5.0, 10.0, 15.0]
results_c = []

println("\n   Варьирование c (β=0.05, γ=0.25):")
for c in cs
    Random.seed!(1234)
    p = [0.05, c, 0.25]
    m = MakeSIRModel(uθ, p)
    activate(m)
    sir_run(m, tmax)
    data = out(m)

    peak_I = maximum(data.I)
    peak_time = data.t[argmax(data.I)]
    final_R = data.R[end]

    push!(results_c, (c=c, peak_I=peak_I, peak_time=peak_time, final_R=final_R))
    println("   c=$c: пик I=$(round(peak_I)), время=$(round(peak_time, digits=2)), переболело=$(round(final_R))")
end

gammas = [0.15, 0.25, 0.35]
results_gamma = []

println("\n   Варьирование γ (β=0.05, c=10.0):")
for γ in gammas
    Random.seed!(1234)
    p = [0.05, 10.0, γ]
    m = MakeSIRModel(uθ, p)
    activate(m)
    sir_run(m, tmax)
    data = out(m)

    peak_I = maximum(data.I)
    peak_time = data.t[argmax(data.I)]
    final_R = data.R[end]

    push!(results_gamma, (γ=γ, peak_I=peak_I, peak_time=peak_time, final_R=final_R))
    println("   γ=$γ: пик I=$(round(peak_I)), время=$(round(peak_time, digits=2)), переболело=$(round(final_R))")
end

p1 = plot(title="Чувствительность к β",
          xlabel="β",
          ylabel="Пик I")
p1 = scatter!([r.β for r in results_beta],
              [r.peak_I for r in results_beta],
              label="Пик I",
              markersize=8)

p2 = plot(title="Чувствительность к c",
          xlabel="c",
          ylabel="Пик I")
p2 = scatter!([r.c for r in results_c],
              [r.peak_I for r in results_c],
              label="Пик I",
              markersize=8)

p3 = plot(title="Чувствительность к γ",
          xlabel="γ",
          ylabel="Пик I")
p3 = scatter!([r.γ for r in results_gamma],
              [r.peak_I for r in results_gamma],
              label="Пик I",
              markersize=8)

plot(p1, p2, p3, layout=(1,3), size=(900, 400))

savefig(plotsdir("sensitivity_analysis.png"))
println("   График чувствительности сохранён: plots/sensitivity_analysis.png")

println("\n3. Сравнение стохастического и детерминированного выздоровления...")

Random.seed!(1234)
m_stoch = MakeSIRModel(uθ, p_default)
activate(m_stoch)
sir_run(m_stoch, tmax)
data_stoch = out(m_stoch)

Random.seed!(1234)
m_det = MakeSIRModel_deterministic_recovery(uθ, p_default)
activate_deterministic(m_det)
sir_run(m_det, tmax)
data_det = out(m_det)

plot(data_stoch.t, data_stoch.I,
     label="Стохастическое выздоровление",
     xlabel="Время",
     ylabel="Инфицированные",
     linewidth=2,
     linestyle=:solid)

plot!(data_det.t, data_det.I,
      label="Детерминированное выздоровление",
      linewidth=2,
      linestyle=:dash)

title!("Сравнение типов выздоровления")

savefig(plotsdir("recovery_comparison.png"))
println("   Сравнительный график сохранён: plots/recovery_comparison.png")

println("\n4. Оценка производительности...")

uθ_large = [9990, 10, 0]
tmax_short = 10.0

println("   Создание модели с N=10000 индивидов...")
@time begin
    m_large = MakeSIRModel(uθ_large, p_default)
    activate(m_large)
end

println("\n   Бенчмарк sir_run для 10000 индивидов:")
bench_result = @benchmark sir_run($m_large, $tmax_short) samples=3

println("   Среднее время: $(round(median(bench_result).time/1e9, digits=3)) секунд")
println("   Минимальное время: $(round(minimum(bench_result).time/1e9, digits=3)) секунд")
println("   Максимальное время: $(round(maximum(bench_result).time/1e9, digits=3)) секунд")

println("\n5. Сохранение результатов в CSV...")

timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")

filename = "sir_$(uθ[1])_$(uθ[2])_$(p_default[1])_$(p_default[2])_$(p_default[3])_$timestamp.csv"

csv_path = datadir("sims", filename)

mkpath(datadir("sims"))

CSV.write(csv_path, data_des)
println("   Результаты сохранены: $csv_path")

println("\n" * "="^60)
println("ИТОГОВАЯ СТАТИСТИКА")
println("="^60)

println("Начальные условия: S₀=$(uθ[1]), I₀=$(uθ[2]), R₀=$(uθ[3])")
println("Параметры: β=$(p_default[1]), c=$(p_default[2]), γ=$(p_default[3])")

R0 = p_default[1] * p_default[2] / p_default[3]
println("R₀ = β·c/γ = $(round(R0, digits=3))")

println("Конечная численность: S=$(data_des.S[end]), I=$(data_des.I[end]), R=$(data_des.R[end])")

println("Максимальное число инфицированных: $(maximum(data_des.I))")
println("Время пика: $(round(data_des.t[argmax(data_des.I)], digits=2))")

final_fraction = data_des.R[end] / sum(uθ) * 100
println("Общая доля переболевших: $(round(final_fraction, digits=1))%")

println("\nРабота завершена успешно!")

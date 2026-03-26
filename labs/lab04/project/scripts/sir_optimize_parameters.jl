using DrWatson
@quickactivate "project"
using BlackBoxOptim, Random, Statistics, JLD2  # ← добавить JLD2
include(srcdir("sir_model.jl"))

# Целевая функция: минимизируем пиковую заболеваемость и смертность
function cost_multi(x)
    # x[1]: β_und, x[2]: detection_time, x[3]: death_rate
    replicates = 5
    peak_vals = Float64[]
    dead_vals = Int[]
    
    for rep in 1:replicates
        model = initialize_sir(;
            Ns = [1000, 1000, 1000],
            β_und = fill(x[1], 3),
            β_det = fill(x[1]/10, 3),
            infection_period = 14,
            detection_time = round(Int, x[2]),  # x[2] - время выявления
            death_rate = x[3],                  # x[3] - смертность
            reinfection_probability = 0.1,
            Is = [0, 0, 1],
            seed = 42 + rep,
            # n_steps НЕ передаем! он не нужен в initialize_sir
        )
        
        infected_frac(model) = count(a.status == :I for a in allagents(model)) / nagents(model)
        peak_infected = 0.0
        
        for step in 1:100
            Agents.step!(model, 1)
            frac = infected_frac(model)
            if frac > peak_infected
                peak_infected = frac
            end
        end
        
        push!(peak_vals, peak_infected)
        push!(dead_vals, 3000 - nagents(model))
    end
    
    return (mean(peak_vals), mean(dead_vals) / 3000)  # доля умерших
end

# Запуск оптимизации
result = bboptimize(
    cost_multi,
    Method = :borg_moea,
    FitnessScheme = ParetoFitnessScheme{2}(is_minimizing = true),
    SearchRange = [
        (0.2, 1.0),      # β_und
        (3.0, 14.0),     # detection_time
        (0.01, 0.1),     # death_rate
    ],
    NumDimensions = 3,
    MaxTime = 120,       # 2 минуты
    TraceMode = :compact,
)

best = best_candidate(result)
fitness = best_fitness(result)

println("Оптимальные параметры:")
println("β_und = $(best[1])")
println("Время выявления = $(round(Int, best[2])) дней")
println("Смертность = $(best[3])")
println("Достигнутые показатели:")
println("Пик заболеваемости: $(fitness[1])")
println("Доля умерших: $(fitness[2])")

# Исправлено: @save, а не save
@save datadir("optimization_result.jld2") best fitness

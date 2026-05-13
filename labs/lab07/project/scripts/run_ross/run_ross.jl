using DrWatson
using Distributions
using ConcurrentSim
using ResumableFunctions
using Random
using StableRNGs
using DataFrames
using Plots
using CSV
using Statistics
using Printf          # ← ДОБАВЛЕНО для @sprintf

const RUNS = 20                    # количество прогонов для статистики
const N_DEFAULT = 10               # работающие машины
const S_DEFAULT = 3                # резервные машины
const R_DEFAULT = 2                # количество ремонтников
const LAMBDA = 100.0               # среднее время безотказной работы (часов)
const MU = 1.0                     # среднее время ремонта (часов)
const SEED = 150

rng = StableRNG(SEED)
failure_dist = Exponential(LAMBDA)
repair_dist = Exponential(MU)

mutable struct Monitor
    time::Vector{Float64}
    operational::Vector{Int}       # работающие машины
    spare::Vector{Int}             # резервные машины
    in_repair::Vector{Int}         # машины в ремонте
    queue_length::Vector{Int}      # очередь на ремонт
end

Monitor() = Monitor(Float64[], Int[], Int[], Int[], Int[])

function record!(mon::Monitor, env, operational, spare, in_repair, queue_len)
    push!(mon.time, now(env))
    push!(mon.operational, operational)
    push!(mon.spare, spare)
    push!(mon.in_repair, in_repair)
    push!(mon.queue_length, queue_len)
end

mutable struct SystemState
    operational::Int
    spare::Int
    in_repair::Int
end

mutable struct RepairStats
    busy::Int      # количество занятых ремонтников
    queue::Int     # длина очереди на ремонт
end

@resumable function machine(
    env::Environment,
    repair_facility::Resource,
    state::SystemState,
    repair_stats::RepairStats,
    mon::Monitor,
    N::Int,
    S::Int,
    R::Int
)
    while true

        @yield timeout(env, rand(rng, failure_dist))

        state.operational -= 1
        state.spare -= 1

        record!(mon, env, state.operational, state.spare,
                state.in_repair, repair_stats.queue)

        if state.spare < 0
            throw(StopSimulation("System crash at time $(now(env)): no spares"))
        end

        if repair_stats.busy >= R
            repair_stats.queue += 1
            record!(mon, env, state.operational, state.spare,
                    state.in_repair, repair_stats.queue)
        end

        @yield request(repair_facility)
        repair_stats.busy += 1
        if repair_stats.queue > 0
            repair_stats.queue -= 1
        end
        state.in_repair += 1

        record!(mon, env, state.operational, state.spare,
                state.in_repair, repair_stats.queue)

        @yield timeout(env, rand(rng, repair_dist))

        @yield release(repair_facility)
        repair_stats.busy -= 1
        state.in_repair -= 1
        state.spare += 1

        if state.spare > 0 && state.operational < N
            state.operational += 1
            state.spare -= 1
        end

        record!(mon, env, state.operational, state.spare,
                state.in_repair, repair_stats.queue)
    end
end

function setup_system(env, repair_facility, state, repair_stats, mon, N, S, R)
    for i in 1:N
        @process machine(env, repair_facility, state, repair_stats, mon, N, S, R)
    end
    return nothing
end

function run_single(N::Int, S::Int, R::Int)
    sim = Simulation()
    repair_facility = Resource(sim, R)
    state = SystemState(N, S, 0)
    repair_stats = RepairStats(0, 0)
    mon = Monitor()

    record!(mon, sim, state.operational, state.spare, state.in_repair, 0)

    setup_system(sim, repair_facility, state, repair_stats, mon, N, S, R)

    msg = run(sim)
    crash_time = now(sim)

    return crash_time, msg, mon
end

function run_experiments()
    results = DataFrame(
        N=Int[], S=Int[], R=Int[],
        mean_time=Float64[], std_time=Float64[],
        min_time=Float64[], max_time=Float64[]
    )

    n_values = [5, 10, 15]
    s_values = [2, 3, 5]
    r_values = [1, 2, 3]

    total_combinations = length(n_values) * length(s_values) * length(r_values)
    current = 0

    for n in n_values
        for s in s_values
            for r in r_values
                current += 1
                println("[$current/$total_combinations] Running N=$n, S=$s, R=$r...")

                times = Float64[]

                for run_idx in 1:RUNS
                    try
                        t, _, _ = run_single(n, s, r)
                        if t > 0
                            push!(times, t)
                        end
                    catch e
                        push!(times, 0.0)
                    end
                end

                times = filter(x -> x > 0, times)

                if length(times) > 0
                    push!(results, (
                        n, s, r,
                        mean(times),
                        std(times),
                        minimum(times),
                        maximum(times)
                    ))
                else
                    push!(results, (n, s, r, 0.0, 0.0, 0.0, 0.0))
                end
            end
        end
    end

    return results
end

function analytical_mttf(N::Int, S::Int, R::Int, lambda::Float64, mu::Float64)

    λ = N / lambda      # суммарная интенсивность отказов
    μ = R / mu          # суммарная интенсивность ремонта

    if λ >= μ
        return Inf      # система нестационарна
    end

    ρ = λ / μ           # загрузка системы

    if R == 1

        return (S + 1) / λ * (1 / (1 - ρ))
    else

        return (S + 1) / λ * (1 / (1 - ρ)) * (1 - ρ^R) / (1 - ρ^(R+1))
    end
end

function plot_results(mon::Monitor, crash_time::Float64, n::Int, s::Int, r::Int)
    if length(mon.time) == 0
        @warn "Нет данных для построения графиков"
        return nothing
    end

    total_good = [mon.operational[i] + mon.spare[i] for i in 1:length(mon.time)]

    p1 = plot(mon.time, total_good,
        label="Исправные машины (всего)",
        xlabel="Время (часы)", ylabel="Количество",
        title="Динамика числа исправных машин (N=$n, S=$s, R=$r)",
        linewidth=2, color=:blue)
    vline!([crash_time], label="Падение системы", linestyle=:dash, linewidth=2, color=:red)
    hline!([n], label="Необходимо для работы", linestyle=:dot, color=:green)

    p2 = plot(mon.time, mon.operational,
        label="Работающие машины",
        xlabel="Время (часы)", ylabel="Количество",
        title="Работающие машины",
        linewidth=2, color=:blue)
    hline!([n], label="Требуется N=$n", linestyle=:dot, color=:red)

    p3 = plot(mon.time, mon.spare,
        label="Резервные машины",
        xlabel="Время (часы)", ylabel="Количество",
        title="Резерв",
        linewidth=2, color=:green)

    p4 = plot(mon.time, mon.queue_length,
        label="Длина очереди на ремонт",
        xlabel="Время (часы)", ylabel="Очередь",
        title="Очередь на ремонт (R=$r ремонтников)",
        linewidth=2, color=:orange)

    p = plot(p1, p2, p3, p4, layout=(4,1), size=(800, 1000))

    savefig("ross_results.png")
    println("Графики сохранены в ross_results.png")

    return p
end

function plot_heatmap(results::DataFrame)

    r1_data = filter(:R => ==(1), results)
    if nrow(r1_data) > 0 && length(unique(r1_data.S)) > 1 && length(unique(r1_data.N)) > 1
        p1 = heatmap(
            unique(r1_data.S),
            unique(r1_data.N),
            reshape(r1_data.mean_time, length(unique(r1_data.N)), length(unique(r1_data.S))),
            xlabel="Резерв (S)", ylabel="Работающие машины (N)",
            title="Среднее время до падения (R=1)",
            color=:viridis,
            aspect_ratio=:equal
        )
        savefig("ross_heatmap_R1.png")
        println("Тепловая карта для R=1 сохранена в ross_heatmap_R1.png")
    end

    r2_data = filter(:R => ==(2), results)
    if nrow(r2_data) > 0 && length(unique(r2_data.S)) > 1 && length(unique(r2_data.N)) > 1
        p2 = heatmap(
            unique(r2_data.S),
            unique(r2_data.N),
            reshape(r2_data.mean_time, length(unique(r2_data.N)), length(unique(r2_data.S))),
            xlabel="Резерв (S)", ylabel="Работающие машины (N)",
            title="Среднее время до падения (R=2)",
            color=:viridis,
            aspect_ratio=:equal
        )
        savefig("ross_heatmap_R2.png")
        println("Тепловая карта для R=2 сохранена в ross_heatmap_R2.png")
    end
end

function compare_with_analytical()
    println("\n" * "="^60)
    println("СРАВНЕНИЕ СИМУЛЯЦИИ С АНАЛИТИКОЙ")
    println("="^60)

    test_cases = [
        (10, 2, 1), (10, 3, 1), (10, 5, 1),
        (10, 2, 2), (10, 3, 2), (10, 5, 2),
        (10, 2, 3), (10, 3, 3), (10, 5, 3)
    ]

    println("\nN       S       R       Симуляция    Аналитика    Отклонение   Доверительный")
    println("-"^80)

    for (n, s, r) in test_cases

        sim_times = Float64[]
        for run_idx in 1:RUNS
            try
                t, _, _ = run_single(n, s, r)
                if t > 0
                    push!(sim_times, t)
                end
            catch

            end
        end

        if length(sim_times) > 0
            sim_mean = mean(sim_times)
            sim_std = std(sim_times)
            sim_se = sim_std / sqrt(length(sim_times))
        else
            sim_mean = 0.0
            sim_std = 0.0
            sim_se = 0.0
        end

        analytic = analytical_mttf(n, s, r, LAMBDA, MU)

        if analytic > 0 && analytic != Inf && sim_mean > 0
            deviation = (sim_mean - analytic) / analytic * 100
            dev_str = @sprintf("%+.1f%%", deviation)
        else
            dev_str = "N/A"
        end

        if sim_se > 0
            ci = @sprintf("%.1f ± %.1f", sim_mean, 1.96 * sim_se)
        else
            ci = "N/A"
        end

        println(@sprintf("%-8d %-8d %-8d %10.1f %10.1f %12s %16s",
                n, s, r, sim_mean, analytic, dev_str, ci))
    end
end

function analyze_repair_utilization()
    println("\n" * "="^60)
    println("АНАЛИЗ ЗАГРУЗКИ РЕМОНТНИКОВ")
    println("="^60)

    configs = [(10, 3, 1), (10, 3, 2), (10, 3, 3)]

    for (n, s, r) in configs
        println("\nКонфигурация: N=$n, S=$s, R=$r")
        println("-"^40)

        try
            _, _, mon = run_single(n, s, r)

            if length(mon.queue_length) > 0
                avg_queue = mean(mon.queue_length)
                max_queue = maximum(mon.queue_length)

                println("  Средняя длина очереди: $(round(avg_queue, digits=2))")
                println("  Максимальная длина очереди: $max_queue")
            else
                println("  Нет данных об очереди")
            end
        catch e
            println("  Ошибка: $e")
        end
    end
end

function main()
    println("="^60)
    println("МОДЕЛЬ РОССА")
    println("="^60)
    println("\nПараметры по умолчанию:")
    println("  N = $N_DEFAULT (работающие машины)")
    println("  S = $S_DEFAULT (резервные машины)")
    println("  R = $R_DEFAULT (ремонтников)")
    println("  λ = 1/$LAMBDA отказов/час (на одну машину)")
    println("  μ = 1/$MU ремонт/час (на один ремонт)")

    total_failure_rate = N_DEFAULT / LAMBDA
    total_repair_rate = R_DEFAULT / MU
    println("\nСуммарные интенсивности:")
    println("  Отказы: $(round(total_failure_rate, digits=3)) отказов/час")
    println("  Ремонт: $(round(total_repair_rate, digits=3)) ремонтов/час")

    if total_failure_rate >= total_repair_rate
        println("  ⚠️  ВНИМАНИЕ: Система нестационарна (отказы > ремонта)")
    else
        println("  ✓ Система стационарна")
    end

    println("\n" * "="^60)
    println("ОДИНОЧНЫЙ ПРОГОН (визуализация)")
    println("="^60)
    println("Запуск симуляции с параметрами по умолчанию...")

    try
        crash_time, msg, mon = run_single(N_DEFAULT, S_DEFAULT, R_DEFAULT)
        println(msg)
        println("Время до падения: $(round(crash_time, digits=2)) часов")

        plot_results(mon, crash_time, N_DEFAULT, S_DEFAULT, R_DEFAULT)
    catch e
        println("Ошибка при запуске: $e")
    end

    println("\n" * "="^60)
    println("МАССОВЫЕ ПРОГОНЫ (RUNS=$RUNS)")
    println("="^60)

    results = run_experiments()

    println("\nРЕЗУЛЬТАТЫ ЭКСПЕРИМЕНТОВ:")
    println("-"^80)
    show(results, allcols=true, eltypes=false)

    if !isdir("results")
        mkdir("results")
    end

    csv_path = joinpath("results", "ross_results.csv")
    CSV.write(csv_path, results)
    println("\n\nРезультаты сохранены в $csv_path")

    txt_path = joinpath("results", "ross_results.txt")
    open(txt_path, "w") do f
        write(f, "РЕЗУЛЬТАТЫ ЭКСПЕРИМЕНТОВ МОДЕЛИ РОССА\n")
        write(f, "="^60 * "\n")
        write(f, "Параметры: RUNS=$RUNS\n\n")
        for row in eachrow(results)
            write(f, @sprintf("N=%2d, S=%2d, R=%2d: %8.2f ± %8.2f часов (min=%8.2f, max=%8.2f)\n",
                row.N, row.S, row.R, row.mean_time, row.std_time, row.min_time, row.max_time))
        end
    end
    println("Текстовые результаты сохранены в $txt_path")

    println("\n" * "="^60)
    println("ВИЗУАЛИЗАЦИЯ РЕЗУЛЬТАТОВ")
    println("="^60)
    plot_heatmap(results)

    compare_with_analytical()

    analyze_repair_utilization()

    println("\n" * "="^60)
    println("АНАЛИЗ ЗАВЕРШЁН")
    println("="^60)

    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

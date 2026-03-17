# Daisyworld: серия экспериментов при изменении солнечной светимости
#
# В данном эксперименте проводится серия симуляций модели Daisyworld
# с различными параметрами. Основная цель — исследовать влияние
# начального количества ромашек и максимального возраста агентов
# на динамику системы при изменяющейся солнечной светимости.

# ## Подключение пакетов

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using DrWatson
@quickactivate "project"

using Agents
using DataFrames
using Plots
using CairoMakie

# ## Подключение реализации модели

include(srcdir("daisyworld.jl"))

# ## Определение типов ромашек

black(a) = a.breed == :black
white(a) = a.breed == :white

# Будем считать количество агентов каждого типа
adata = [(black, count), (white, count)]

# ## Параметры эксперимента
#
# Некоторые параметры задаются как массивы значений.
# Для них будут автоматически сформированы комбинации
# параметров, которые будут протестированы.

param_dict = Dict(
    :griddims => (30, 30),
    :max_age => [25, 40],
    :init_white => [0.2, 0.8],
    :init_black => 0.2,
    :albedo_white => 0.75,
    :albedo_black => 0.25,
    :surface_albedo => 0.4,
    :solar_change => 0.005,
    :solar_luminosity => 1.0,
    :scenario => :ramp,
    :seed => 165,
)

# ## Формирование списка параметров

params_list = dict_list(param_dict)

# ## Проведение серии симуляций

for params in params_list

    model = daisyworld(; params...)

    # Средняя температура поверхности
    temperature(model) = StatsBase.mean(model.temperature)

    # Сохраняем температуру и солнечную светимость
    mdata = [temperature, :solar_luminosity]

    agent_df, model_df = run!(
        model,
        1000;
        adata = adata,
        mdata = mdata
    )

    # ## Построение графиков

    figure = CairoMakie.Figure(size = (600, 600))

    # --- динамика численности ромашек ---

    ax1 = figure[1, 1] = Axis(
        figure,
        ylabel = "daisy count"
    )

    black1 = lines!(
        ax1,
        agent_df[!, :time],
        agent_df[!, :count_black],
        color = :red
    )

    white1 = lines!(
        ax1,
        agent_df[!, :time],
        agent_df[!, :count_white],
        color = :blue
    )

    figure[1, 2] = Legend(
        figure,
        [black1, white1],
        ["black", "white"]
    )

    # --- температура планеты ---

    ax2 = figure[2, 1] = Axis(
        figure,
        ylabel = "temperature"
    )

    lines!(
        ax2,
        model_df[!, :time],
        model_df[!, :temperature],
        color = :red
    )

    # --- солнечная светимость ---

    ax3 = figure[3, 1] = Axis(
        figure,
        xlabel = "tick",
        ylabel = "luminosity"
    )

    lines!(
        ax3,
        model_df[!, :time],
        model_df[!, :solar_luminosity],
        color = :red
    )

    # скрываем подписи оси X у верхних графиков
    for ax in (ax1, ax2)
        ax.xticklabelsvisible = false
    end

    # ## Сохранение результата

    plt_name = savename("daisy-luminosity", params) * ".png"

    save(plotsdir(plt_name), figure)
end


# ## Анализ результатов
#
# Для каждой комбинации параметров выполняется симуляция модели
# Daisyworld и строятся три графика:
#
# 1. изменение количества черных и белых ромашек
# 2. изменение средней температуры поверхности
# 3. изменение солнечной светимости
#
# Поскольку используется сценарий :ramp, солнечная светимость
# постепенно изменяется со временем. Это приводит к изменению
# температуры среды, что влияет на динамику популяций ромашек.
#
# Сравнение графиков для различных параметров позволяет оценить,
# как максимальный возраст ромашек и их начальное распределение
# влияют на устойчивость экосистемы и её способность
# адаптироваться к изменению внешних условий.

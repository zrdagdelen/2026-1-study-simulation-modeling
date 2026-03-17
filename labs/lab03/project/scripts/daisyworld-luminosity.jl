# Daisyworld: влияние солнечной светимости
#
# В этом эксперименте исследуется влияние изменения
# солнечной светимости на динамику популяции ромашек
# и температуру планеты в модели Daisyworld.

# ## Подключение пакетов и активация проекта

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

# ## Определение типов агентов

black(a) = a.breed == :black
white(a) = a.breed == :white

# Сбор статистики по количеству агентов

adata = [(black, count), (white, count)]

# ## Создание модели
#
# Используется сценарий :ramp — в нём солнечная
# светимость постепенно изменяется во времени.

model = daisyworld(
    solar_luminosity = 1.0,
    scenario = :ramp
)

# ## Сбор данных модели
#
# Будем сохранять:
# - среднюю температуру поверхности
# - текущую солнечную светимость

temperature(model) = StatsBase.mean(model.temperature)

mdata = [temperature, :solar_luminosity]

# ## Запуск симуляции

agent_df, model_df = run!(
    model,
    1000;
    adata = adata,
    mdata = mdata
)

# agent_df содержит данные об агентах
# model_df содержит глобальные параметры модели

# ## Построение графиков

figure = CairoMakie.Figure(size = (600, 600))

# --- График численности ромашек ---

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

# --- График температуры ---

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

# --- График солнечной светимости ---

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

# скрываем подписи оси X на верхних графиках

for ax in (ax1, ax2)
    ax.xticklabelsvisible = false
end

figure

# ## Сохранение графика

save(plotsdir("daisy_luminosity.png"), figure)

# ## Анализ результатов
#
# На графиках представлены:
# 1. изменение количества черных и белых ромашек
# 2. изменение средней температуры поверхности
# 3. изменение солнечной светимости
#
# При постепенном увеличении светимости температура
# планеты также изменяется. Популяции ромашек реагируют
# на эти изменения, поскольку различные типы растений
# имеют разное альбедо и по-разному влияют на нагрев
# поверхности.
#
# Таким образом модель демонстрирует механизм
# климатической саморегуляции: изменение состава
# растительности может компенсировать изменения
# внешних условий.

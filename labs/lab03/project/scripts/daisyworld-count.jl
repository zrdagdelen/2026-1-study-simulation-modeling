# Daisyworld: динамика количества ромашек
#
# В этом скрипте моделируется система Daisyworld и
# анализируется изменение количества черных и белых
# ромашек во времени.

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
#
# В модели есть два типа ромашек:
# - черные (black)
# - белые (white)

black(a) = a.breed == :black
white(a) = a.breed == :white

# Будем считать количество агентов каждого типа
adata = [(black, count), (white, count)]

# ## Создание модели
#
# Параметр solar_luminosity задаёт интенсивность
# солнечного излучения.

model = daisyworld(; solar_luminosity = 1.0)

# ## Запуск симуляции
#
# Функция run! выполняет модель в течение заданного
# числа шагов и собирает статистику.

agent_df, model_df = run!(model, 1000; adata)

# agent_df содержит данные о количестве ромашек
# каждого типа на каждом шаге модели.

# ## Построение графика

figure = Figure(size = (600, 400))

ax = figure[1, 1] = Axis(
    figure,
    xlabel = "tick",
    ylabel = "daisy count"
)

blackl = lines!(
    ax,
    agent_df[!, :time],
    agent_df[!, :count_black],
    color = :black
)

whitel = lines!(
    ax,
    agent_df[!, :time],
    agent_df[!, :count_white],
    color = :orange
)

Legend(
    figure[1, 2],
    [blackl, whitel],
    ["black", "white"],
    labelsize = 12
)

figure

# ## Сохранение результата

save(plotsdir("daisy_count.png"), figure)

# ## Анализ
#
# График показывает изменение количества черных и белых
# ромашек во времени. В зависимости от температуры среды
# один из типов может иметь преимущество в размножении.
#
# Поскольку разные типы ромашек имеют различное альбедо,
# их распространение влияет на температуру поверхности.
# Это приводит к обратной связи между растительностью
# и климатом.
#
# Модель демонстрирует механизм саморегуляции:
# изменение численности растений способствует
# стабилизации температурных условий.

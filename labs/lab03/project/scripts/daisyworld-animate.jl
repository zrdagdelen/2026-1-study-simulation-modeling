# Daisyworld: моделирование и запись видео симуляции
#
# В этом скрипте создаётся агентная модель Daisyworld
# и записывается видео эволюции системы.

# ## Подключение пакетов и активация проекта

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))   # активируем проект

using DrWatson
@quickactivate "project"

using Agents
using DataFrames
using Plots
using CairoMakie

# ## Подключение реализации модели

include(srcdir("daisyworld.jl"))

# ## Создание модели

model = daisyworld()

# В модели каждый агент — ромашка (Daisy).
# Поле breed определяет тип ромашки и используется
# для выбора цвета при визуализации.

daisycolor(a::Daisy) = a.breed

# ## Параметры визуализации

plotkwargs = (
    agent_color = daisycolor,
    agent_size = 20,
    agent_marker = '✿',
    heatarray = :temperature,
    heatkwargs = (colorrange = (-20, 60),),
)

# ## Запуск симуляции и запись видео

abmvideo(
    plotsdir("simulation.mp4"),
    model;
    title = "Daisy World",
    frames = 60,
    plotkwargs...,
)

# ## Анализ
#
# Видео демонстрирует динамику модели Daisyworld.
# Ромашки изменяют альбедо поверхности, что влияет
# на температуру окружающей среды.
#
# В процессе симуляции можно наблюдать:
# - распространение ромашек
# - изменение температурного поля
# - формирование устойчивых областей.
#
# Модель иллюстрирует механизм саморегуляции:
# взаимодействие между растениями и климатом
# может стабилизировать температуру системы.

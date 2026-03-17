using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using DrWatson
@quickactivate "project"  # активирует ваш проект

import Pkg

# Список всех нужных пакетов
packages = [
    "Agents",
    "StatsBase",
    "DataFrames",
    "Plots",
    "CairoMakie",
    "Random"
]

println("Установка пакетов в проект...")
for pkg in packages
    println("Установка $pkg...")
    Pkg.add(pkg)
end

println("Все пакеты установлены!")
Pkg.status()  # покажет установленные версии

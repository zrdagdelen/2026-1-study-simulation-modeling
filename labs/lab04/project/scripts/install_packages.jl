using DrWatson
@quickactivate "project"

using Pkg
Pkg.add.([
    "Distributions",    # для Poisson распределения (нужен в transmit!)
    "BlackBoxOptim",    # для оптимизации
    "CSV",              # для CSV файлов
    "JLD2",             # для JLD2 файлов
    "Literate",         # для литературного кода
    "Quarto"            # для документации
])

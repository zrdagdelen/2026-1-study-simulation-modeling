# ## Итоговый отчёт по модели SIR

# ### Загрузка данных

# Загружаем ранее сохранённые результаты симуляций.

using DrWatson
@quickactivate "project"
using DataFrames, CSV, Plots

# Загрузка данных из CSV-файлов:
# - `sir_det.csv` — детерминированная симуляция
# - `sir_stoch.csv` — стохастическая симуляция
# - `sir_scan.csv` — результаты сканирования параметра β

df_det = CSV.read(datadir("sir_det.csv"), DataFrame)
df_stoch = CSV.read(datadir("sir_stoch.csv"), DataFrame)
df_scan = CSV.read(datadir("sir_scan.csv"), DataFrame)

# ### Сравнение детерминированной и стохастической динамики

# Строим сравнительный график числа инфицированных I(t)
# для детерминированного и стохастического случаев.

n_points = min(length(df_det.time), length(df_stoch.time))
p1 = plot(
    df_det.time[1:n_points],
    [df_det.I[1:n_points] df_stoch.I[1:n_points]],
    label = ["Deterministic I" "Stochastic I"],
    xlabel = "Time",
    ylabel = "Infected",
    title = "Comparison of Deterministic and Stochastic Dynamics",
    linewidth = 2,
)
savefig(plotsdir("comparison.png"))

# ### Анализ чувствительности к параметру β

# Строим график зависимости пика эпидемии (максимального числа инфицированных)
# от коэффициента заражения β.

p2 = plot(
    df_scan.β,
    df_scan.peak_I,
    marker = :circle,
    linewidth = 2,
    xlabel = "β (infection rate)",
    ylabel = "Peak Infected",
    title = "Sensitivity Analysis: Peak I vs β",
)
savefig(plotsdir("sensitivity.png"))

# ### Завершение

println("Отчётные графики сохранены в plots/")

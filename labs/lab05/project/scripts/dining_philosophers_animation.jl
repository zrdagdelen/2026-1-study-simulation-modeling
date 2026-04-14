# # Анимация работы сети Петри для задачи об обедающих философах
#
# ## Введение
#
# Данный скрипт создаёт анимацию, показывающую динамику изменения маркировки
# в классической сети Петри. Анимация позволяет наглядно увидеть, как фишки
# перемещаются между позициями и как возникает deadlock.
#
# ## Загрузка пакетов

using DrWatson
@quickactivate "project"

include(srcdir("DiningPhilosophers.jl"))
using .DiningPhilosophers
using Plots, Random

# ## Параметры
#
# - `N = 3` — количество философов (меньше для лучшей визуализации)
# - `tmax = 30.0` — время моделирования

N = 3
tmax = 30.0

# ## Построение классической сети

net, u0, names = build_classical_network(N)

# ## Стохастическая симуляция
#
# Фиксируем seed для воспроизводимости результатов.

Random.seed!(123)
df = simulate_stochastic(net, u0, tmax)

# ## Создание анимации
#
# Каждый кадр анимации представляет собой столбчатую диаграмму текущей маркировки.
# По оси X — позиции сети, по оси Y — количество фишек.

anim = @animate for row in eachrow(df)
    u = [row[col] for col in propertynames(row) if col != :time]
    bar(
        1:length(u), u,
        legend = false,
        ylims = (0, maximum(u0) + 1),
        xlabel = "Позиция",
        ylabel = "Фишки",
        title = "Время = $(round(row.time, digits=2))",
    )
    xticks!(1:length(u), string.(names), rotation = 45)
end

# ## Сохранение анимации
#
# Анимация сохраняется в формате GIF.

gif(anim, plotsdir("philosophers_simulation.gif"), fps = 2)
println("Анимация сохранена в plots/philosophers_simulation.gif")

# ## Интерпретация
#
# На анимации можно наблюдать:
# - Перемещение фишек между позициями Think, Hungry, Eat и Fork
# - В определённый момент в классической модели наступает deadlock —
#   маркировка перестаёт изменяться, все философы застревают в состоянии Hungry
#
# Это наглядно демонстрирует проблему взаимной блокировки в параллельных системах.

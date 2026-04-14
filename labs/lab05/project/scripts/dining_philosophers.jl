# # Моделирование задачи об обедающих философах с помощью сетей Петри
#
# ## Введение
#
# Данный скрипт выполняет стохастическое моделирование двух вариантов сети Петри
# для классической задачи об обедающих философах:
# 1. Классическая модель (без арбитра) — подвержена deadlock
# 2. Модифицированная модель с арбитром — предотвращает deadlock
#
# ## Загрузка необходимых пакетов

using DrWatson
@quickactivate "project"

include(srcdir("DiningPhilosophers.jl"))
using .DiningPhilosophers
using DataFrames, CSV, Plots

# ## Параметры моделирования
#
# - `N = 5` — количество философов
# - `tmax = 50.0` — время моделирования

N = 5
tmax = 50.0

# ## Моделирование классической сети (без арбитра)
#
# В классической постановке каждый философ пытается взять сначала левую вилку,
# затем правую. Это может привести к взаимной блокировке (deadlock).

println("=== Классическая сеть (без арбитра) ===")
net_classic, u0_classic, _ = build_classical_network(N)

# ### Стохастическая симуляция
#
# Используется алгоритм Гиллеспи для стохастического моделирования.

df_classic = simulate_stochastic(net_classic, u0_classic, tmax)

# ### Сохранение результатов
#
# Результаты сохраняются в формате CSV для дальнейшего анализа.

CSV.write(datadir("dining_classic.csv"), df_classic)

# ### Проверка наличия deadlock

dead = detect_deadlock(df_classic, net_classic)
println("Deadlock обнаружен: $dead")

# ### Визуализация эволюции маркировки
#
# Строятся графики для всех позиций сети: Think, Hungry, Eat, Fork.

plot_classic = plot_marking_evolution(df_classic, N)
savefig(plotsdir("classic_simulation.png"))

# ## Моделирование сети с арбитром
#
# В этой модификации добавлена дополнительная позиция-арбитр с N-1 фишками,
# что предотвращает одновременный захват вилок всеми философами.

println("\n=== Сеть с арбитром ===")
net_arb, u0_arb, _ = build_arbiter_network(N)

# ### Стохастическая симуляция

df_arb = simulate_stochastic(net_arb, u0_arb, tmax)

# ### Сохранение результатов

CSV.write(datadir("dining_arbiter.csv"), df_arb)

# ### Проверка наличия deadlock

dead_arb = detect_deadlock(df_arb, net_arb)
println("Deadlock обнаружен: $dead_arb")

# ### Визуализация эволюции маркировки

plot_arb = plot_marking_evolution(df_arb, N)
savefig(plotsdir("arbiter_simulation.png"))

# ## Выводы
#
# - Классическая сеть демонстрирует deadlock (философы застревают в состоянии Hungry)
# - Сеть с арбитром успешно предотвращает deadlock (система продолжает работать)
#
# Результаты сохранены в:
# - `data/dining_classic.csv` и `data/dining_arbiter.csv` — траектории состояний
# - `plots/classic_simulation.png` и `plots/arbiter_simulation.png` — графики

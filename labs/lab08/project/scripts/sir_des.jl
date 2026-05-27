# # Обзор sir_des
# 
# Данный скрипт является точкой входа для запуска дискретно-событийной
# SIR-модели. Он выполняет:
# 
# 1. Базовый запуск модели с параметрами по умолчанию
# 2. Анализ чувствительности к параметрам `β, c, γ`
# 3. Сравнение стохастической и детерминированной версий выздоровления
# 4. Оценку производительности для популяции из 10000 индивидов
# 5. Сохранение результатов в CSV файл
# 6. Вывод итоговой статистики
# 
# ## Подключение библиотек и модулей

# ### Активация окружения проекта DrWatson
# 
# DrWatson обеспечивает воспроизводимость и управление путями к файлам.

#!/usr/bin/env julia

using DrWatson
@quickactivate "project"

# ### Подключение ядра модели
# 
# Импортируем функции из `sir_model.jl`:
# - `MakeSIRModel` - создание модели
# - `MakeSIRModel_deterministic_recovery` - модель с детерминированным выздоровлением
# - `activate` / `activate_deterministic` - активация процессов
# - `sir_run` - запуск симуляции
# - `out` - сбор результатов

include(srcdir("sir_model.jl"))

# ### Импорт дополнительных библиотек

using Random      # Генерация случайных чисел и фиксация seed
using StatsPlots  # Построение графиков
using BenchmarkTools # Оценка производительности
using CSV         # Сохранение результатов в CSV
using Dates       # Работа с датами для именования файлов

# ## Параметры модели по умолчанию

# Определим основные параметры моделирования:

# **Длительность симуляции** - 40 единиц времени
tmax = 40.0

# **Начальные условия** - 990 восприимчивых, 10 инфицированных, 0 переболевших
uθ = [990, 10, 0]  # S, I, R

# **Параметры модели** - вероятность заражения β=0.05, частота контактов c=10.0, интенсивность выздоровления γ=0.25
p_default = [0.05, 10.0, 0.25]  # β, c, γ

# **Фиксация seed для воспроизводимости**
Random.seed!(1234)

# Вывод заголовка программы
println("="^60)
println("ДИСКРЕТНО-СОБЫТИЙНАЯ SIR МОДЕЛЬ")
println("="^60)

# ## Базовый запуск модели

# Сначала выполним базовый запуск с параметрами по умолчанию.
# Это позволит получить эталонную динамику для сравнения.

println("\n1. Базовый запуск модели...")

# Создание модели
des_model = MakeSIRModel(uθ, p_default)

# Активация процессов агентов
activate(des_model)

# Запуск симуляции до времени tmax
sir_run(des_model, tmax)

# Сбор результатов в DataFrame
data_des = out(des_model)

# ### Визуализация базового запуска
# 
# Построим график изменения численности S, I, R во времени.

plot(data_des.t, [data_des.S data_des.I data_des.R],
     labels=["S" "I" "R"],
     xlabel="Время", 
     ylabel="Численность",
     title="Дискретно-событийная SIR модель (β=0.05, c=10.0, γ=0.25)",
     linewidth=2)

# Сохранение графика в директорию plots/
savefig(plotsdir("sir_des.png"))

println("   График сохранён: plots/sir_des.png")

# ## Анализ чувствительности к параметрам (Задание 8.5.1)
# 
# Проведём анализ влияния каждого параметра на динамику эпидемии.
# Для каждого набора параметров записываем:
# - пик инфицированных
# - время достижения пика
# - итоговую долю переболевших

println("\n2. Анализ чувствительности к параметрам...")

# ### Варьирование $\beta$ (вероятность заражения)
# 
# $\beta$ - вероятность передачи инфекции при контакте.
# Значения варьируются от 0.03 до 0.07.

betas = [0.03, 0.05, 0.07]
results_beta = []

println("\n   Варьирование β (c=10.0, γ=0.25):")
for β in betas
    # Фиксируем seed для воспроизводимости при каждом запуске
    Random.seed!(1234)
    p = [β, 10.0, 0.25]
    m = MakeSIRModel(uθ, p)
    activate(m)
    sir_run(m, tmax)
    data = out(m)
    
    # Вычисление ключевых метрик
    peak_I = maximum(data.I)                        # Пик инфицированных
    peak_time = data.t[argmax(data.I)]              # Время пика
    final_R = data.R[end]                           # Итоговое число переболевших
    
    push!(results_beta, (β=β, peak_I=peak_I, peak_time=peak_time, final_R=final_R))
    println("   β=$β: пик I=$(round(peak_I)), время=$(round(peak_time, digits=2)), переболело=$(round(final_R))")
end

# ### Варьирование c (частота контактов)
# 
# $c$ - среднее число контактов индивида в единицу времени.
# Значения варьируются от 5.0 до 15.0.

cs = [5.0, 10.0, 15.0]
results_c = []

println("\n   Варьирование c (β=0.05, γ=0.25):")
for c in cs
    Random.seed!(1234)
    p = [0.05, c, 0.25]
    m = MakeSIRModel(uθ, p)
    activate(m)
    sir_run(m, tmax)
    data = out(m)
    
    peak_I = maximum(data.I)
    peak_time = data.t[argmax(data.I)]
    final_R = data.R[end]
    
    push!(results_c, (c=c, peak_I=peak_I, peak_time=peak_time, final_R=final_R))
    println("   c=$c: пик I=$(round(peak_I)), время=$(round(peak_time, digits=2)), переболело=$(round(final_R))")
end

# ### Варьирование $\gamma$  (интенсивность выздоровления)
# 
# $\gamma$ - скорость выздоровления (обратная среднему времени болезни).
# Значения варьируются от 0.15 до 0.35.

gammas = [0.15, 0.25, 0.35]
results_gamma = []

println("\n   Варьирование γ (β=0.05, c=10.0):")
for γ in gammas
    Random.seed!(1234)
    p = [0.05, 10.0, γ]
    m = MakeSIRModel(uθ, p)
    activate(m)
    sir_run(m, tmax)
    data = out(m)
    
    peak_I = maximum(data.I)
    peak_time = data.t[argmax(data.I)]
    final_R = data.R[end]
    
    push!(results_gamma, (γ=γ, peak_I=peak_I, peak_time=peak_time, final_R=final_R))
    println("   γ=$γ: пик I=$(round(peak_I)), время=$(round(peak_time, digits=2)), переболело=$(round(final_R))")
end

# ### Визуализация результатов чувствительности
# 
# Построим графики зависимости пика инфицированных от каждого параметра.

# График для $\beta$ 
p1 = plot(title="Чувствительность к β", 
          xlabel="β", 
          ylabel="Пик I")
p1 = scatter!([r.β for r in results_beta], 
              [r.peak_I for r in results_beta], 
              label="Пик I", 
              markersize=8)

# График для $c$
p2 = plot(title="Чувствительность к c", 
          xlabel="c", 
          ylabel="Пик I")
p2 = scatter!([r.c for r in results_c], 
              [r.peak_I for r in results_c], 
              label="Пик I", 
              markersize=8)

# График для $\gamma$ 
p3 = plot(title="Чувствительность к γ", 
          xlabel="γ", 
          ylabel="Пик I")
p3 = scatter!([r.γ for r in results_gamma], 
              [r.peak_I for r in results_gamma], 
              label="Пик I", 
              markersize=8)

# Объединение графиков в один ряд
plot(p1, p2, p3, layout=(1,3), size=(900, 400))

# Сохранение
savefig(plotsdir("sensitivity_analysis.png"))
println("   График чувствительности сохранён: plots/sensitivity_analysis.png")

# ## Сравнение стохастической и детерминированной версий (Задание 8.5.2)
# 
# Сравним два варианта модели:
# - **Стохастическая версия**: время выздоровления ~ Exponential($1/\gamma$)
# - **Детерминированная версия**: время выздоровления = $1/\gamma$ (фиксированное)

println("\n3. Сравнение стохастического и детерминированного выздоровления...")

# ### Стохастическая версия
Random.seed!(1234)
m_stoch = MakeSIRModel(uθ, p_default)
activate(m_stoch)
sir_run(m_stoch, tmax)
data_stoch = out(m_stoch)

# ### Детерминированная версия
Random.seed!(1234)
m_det = MakeSIRModel_deterministic_recovery(uθ, p_default)
activate_deterministic(m_det)
sir_run(m_det, tmax)
data_det = out(m_det)

# ### Сравнительный график
# 
# Наложим кривые инфицированных для обеих версий.

plot(data_stoch.t, data_stoch.I, 
     label="Стохастическое выздоровление", 
     xlabel="Время", 
     ylabel="Инфицированные", 
     linewidth=2, 
     linestyle=:solid)

plot!(data_det.t, data_det.I, 
      label="Детерминированное выздоровление", 
      linewidth=2, 
      linestyle=:dash)

title!("Сравнение типов выздоровления")

savefig(plotsdir("recovery_comparison.png"))
println("   Сравнительный график сохранён: plots/recovery_comparison.png")

# ## Оценка производительности (Задание 8.5.3)
# 
# Измерим время выполнения модели для популяции из 10000 индивидов.
# Используем макрос `@benchmark` из пакета BenchmarkTools.

println("\n4. Оценка производительности...")

# ### Создание модели большого размера
# 
# Популяция: 9990 восприимчивых + 10 инфицированных = 10000 индивидов
# Время симуляции сокращено до 10 единиц для ускорения

uθ_large = [9990, 10, 0]
tmax_short = 10.0

println("   Создание модели с N=10000 индивидов...")
@time begin
    m_large = MakeSIRModel(uθ_large, p_default)
    activate(m_large)
end

# ### Бенчмаркинг
# 
# Выполняем 3 прогона для статистической оценки

println("\n   Бенчмарк sir_run для 10000 индивидов:")
bench_result = @benchmark sir_run($m_large, $tmax_short) samples=3

# Вывод результатов в секундах
println("   Среднее время: $(round(median(bench_result).time/1e9, digits=3)) секунд")
println("   Минимальное время: $(round(minimum(bench_result).time/1e9, digits=3)) секунд")
println("   Максимальное время: $(round(maximum(bench_result).time/1e9, digits=3)) секунд")

# ## Сохранение результатов в CSV (Задание 8.5.4)
# 
# Автоматическое сохранение итоговой таблицы с уникальным именем,
# содержащим параметры запуска.

println("\n5. Сохранение результатов в CSV...")

# Генерация временной метки
timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")

# Формирование имени файла с параметрами
filename = "sir_$(uθ[1])_$(uθ[2])_$(p_default[1])_$(p_default[2])_$(p_default[3])_$timestamp.csv"

# Путь для сохранения
csv_path = datadir("sims", filename)

# Создание директории (если не существует)
mkpath(datadir("sims"))

# Сохранение DataFrame в CSV
CSV.write(csv_path, data_des)
println("   Результаты сохранены: $csv_path")

# ## Итоговая статистика
# 
# Вывод сводных результатов моделирования.

println("\n" * "="^60)
println("ИТОГОВАЯ СТАТИСТИКА")
println("="^60)

# Основные параметры
println("Начальные условия: S₀=$(uθ[1]), I₀=$(uθ[2]), R₀=$(uθ[3])")
println("Параметры: β=$(p_default[1]), c=$(p_default[2]), γ=$(p_default[3])")

# Репродуктивное число
R0 = p_default[1] * p_default[2] / p_default[3]
println("R₀ = β·c/γ = $(round(R0, digits=3))")

# Конечные значения
println("Конечная численность: S=$(data_des.S[end]), I=$(data_des.I[end]), R=$(data_des.R[end])")

# Пик эпидемии
println("Максимальное число инфицированных: $(maximum(data_des.I))")
println("Время пика: $(round(data_des.t[argmax(data_des.I)], digits=2))")

# Итоговая доля переболевших
final_fraction = data_des.R[end] / sum(uθ) * 100
println("Общая доля переболевших: $(round(final_fraction, digits=1))%")

println("\nРабота завершена успешно!")

# ## Выводы
# 
# На основе выполненного анализа можно сделать следующие выводы:
# 
# 1. **Чувствительность к параметрам**:
#    - Увеличение  $\beta$ и c приводит к росту пика инфицированных и уменьшению времени пика
#    - Увеличение  $\gamma$ снижает пик и замедляет распространение эпидемии
# 
# 2. **Сравнение версий**:
#    - Детерминированная версия дает более гладкую кривую без флуктуаций
#    - Стохастическая версия более реалистична для малых популяций
# 
# 3. **Производительность**:
#    - Для популяции 10000 индивидов модель выполняет симуляцию за приблизительно 0.02 секунды (время создания модели)
# 
# 4. **Репродуктивное число**:
#    - При $R_0$ > 1 эпидемия распространяется
#    - При $R_0$ < 1 эпидемия затухает
# 

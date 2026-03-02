# # Исследование влияния параметров модели SIR
#
# ## Введение
# В этом исследовании мы анализируем, как различные параметры модели SIR
# влияют на динамику эпидемии. Модель описывается системой уравнений:
#
# $$ \frac{dS}{dt} = -\beta \frac{IS}{N} $$
# $$ \frac{dI}{dt} = \beta \frac{IS}{N} - \gamma I $$
# $$ \frac{dR}{dt} = \gamma I $$
#
# ## Подключение пакетов

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using DifferentialEquations
using DataFrames
using StatsPlots
using Plots
using LaTeXStrings

# ## Определение функции модели

function sir_ode!(du, u, p, t)
    (S, I, R) = u
    (β, c, γ) = p
    N = S + I + R
    du[1] = -β * c * I / N * S
    du[2] = β * c * I / N * S - γ * I
    du[3] = γ * I
end

# Начальные условия (одинаковые для всех экспериментов)
u0 = [990.0, 10.0, 0.0]
tspan = (0.0, 40.0)

# ## Вспомогательная функция для запуска модели

function run_sir(p; title="")
    prob = ODEProblem(sir_ode!, u0, tspan, p)
    sol = solve(prob, dt=0.1)
    
    df = DataFrame(Tables.table(sol'))
    rename!(df, ["S", "I", "R"])
    df[!, :t] = sol.t
    df[!, :N] = df.S + df.I + df.R
    
    R0 = (p[2] * p[1]) / p[3]
    return df, R0
end

# # Эксперимент 1: Влияние вероятности заражения (β)
#
# Параметр β отвечает за вероятность передачи инфекции при контакте.
# Чем выше β, тем заразнее болезнь.

println("\n" * "="^60)
println("ИССЛЕДОВАНИЕ 1: Влияние вероятности заражения (β)")
println("="^60)

β_values = [0.02, 0.05, 0.10, 0.15]
c_fixed = 10.0
γ_fixed = 0.25

plt_β = plot(title="Влияние вероятности заражения β", 
             xlabel="Время, дни", 
             ylabel="Инфицированные", 
             legend=:topright)

for β in β_values
    p = [β, c_fixed, γ_fixed]
    df, R0 = run_sir(p)
    plot!(plt_β, df.t, df.I, label="β=$(β), R0=$(round(R0, digits=2))", linewidth=2)
    
    peak_idx = argmax(df.I)
    println("β = $(β): R0 = $(round(R0, digits=2)), пик = $(round(df.I[peak_idx], digits=1)) на $(round(df.t[peak_idx], digits=1)) день")
end
plt_β

# **Анализ:** С увеличением β растет R0, эпидемия развивается быстрее,
# пик становится выше и наступает раньше.

# # Эксперимент 2: Влияние числа контактов (c)
#
# Параметр c отражает социальную активность - среднее число контактов
# одного человека в день. Это поведенческий фактор.

println("\n" * "="^60)
println("ИССЛЕДОВАНИЕ 2: Влияние числа контактов (c)")
println("="^60)

c_values = [5.0, 10.0, 15.0, 20.0]
β_fixed = 0.05
γ_fixed = 0.25

plt_c = plot(title="Влияние числа контактов c", 
             xlabel="Время, дни", 
             ylabel="Инфицированные", 
             legend=:topright)

for c in c_values
    p = [β_fixed, c, γ_fixed]
    df, R0 = run_sir(p)
    plot!(plt_c, df.t, df.I, label="c=$(c), R0=$(round(R0, digits=2))", linewidth=2)
    
    peak_idx = argmax(df.I)
    println("c = $(c): R0 = $(round(R0, digits=2)), пик = $(round(df.I[peak_idx], digits=1)) на $(round(df.t[peak_idx], digits=1)) день")
end
plt_c

# **Анализ:** Число контактов линейно влияет на R0. Снижение контактов
# (карантин, самоизоляция) эффективно замедляет эпидемию.

# # Эксперимент 3: Влияние скорости выздоровления (γ)
#
# Параметр γ определяет, как быстро люди выздоравливают. 
# Средняя продолжительность болезни = 1/γ.

println("\n" * "="^60)
println("ИССЛЕДОВАНИЕ 3: Влияние скорости выздоровления (γ)")
println("="^60)

γ_values = [0.1, 0.2, 0.25, 0.33]  # 1/γ = 10, 5, 4, 3 дней
β_fixed = 0.05
c_fixed = 10.0

plt_γ = plot(title="Влияние скорости выздоровления γ", 
             xlabel="Время, дни", 
             ylabel="Инфицированные", 
             legend=:topright)

for γ in γ_values
    p = [β_fixed, c_fixed, γ]
    df, R0 = run_sir(p)
    plot!(plt_γ, df.t, df.I, label="γ=$(γ) (болезнь $(round(1/γ, digits=1)) дней), R0=$(round(R0, digits=2))", linewidth=2)
    
    peak_idx = argmax(df.I)
    println("γ = $(γ) (1/γ = $(round(1/γ, digits=1)) дней): R0 = $(round(R0, digits=2)), пик = $(round(df.I[peak_idx], digits=1)) на $(round(df.t[peak_idx], digits=1)) день")
end
plt_γ

# **Анализ:** Увеличение γ (более быстрое выздоровление) снижает R0.
# Лечение и изоляция больных уменьшают время заразности.

# # Эксперимент 4: Сравнение мер контроля
#
# Сравним эффективность разных стратегий борьбы с эпидемией:
# - Карантин (снижение контактов)
# - Маски и гигиена (снижение вероятности заражения)
# - Быстрое лечение (ускорение выздоровления)

println("\n" * "="^60)
println("ИССЛЕДОВАНИЕ 4: Сравнение мер контроля")
println("="^60)

# Базовый сценарий (без мер)
p_base = [0.05, 10.0, 0.25]
df_base, R0_base = run_sir(p_base)

# Сценарий 1: Снижение контактов (карантин) - c уменьшено на 50%
p_quarantine = [0.05, 5.0, 0.25]
df_quarantine, R0_quarantine = run_sir(p_quarantine)

# Сценарий 2: Защита (маски) - β уменьшено на 50%
p_masks = [0.025, 10.0, 0.25]
df_masks, R0_masks = run_sir(p_masks)

# Сценарий 3: Быстрое лечение - γ увеличено на 50%
p_treatment = [0.05, 10.0, 0.375]
df_treatment, R0_treatment = run_sir(p_treatment)

plt_control = plot(title="Сравнение мер контроля эпидемии", 
                   xlabel="Время, дни", 
                   ylabel="Инфицированные", 
                   legend=:topright, 
                   linewidth=2)

plot!(plt_control, df_base.t, df_base.I, label="Без мер: R0=$(round(R0_base, digits=2))")
plot!(plt_control, df_quarantine.t, df_quarantine.I, label="Карантин (c↓50%): R0=$(round(R0_quarantine, digits=2))")
plot!(plt_control, df_masks.t, df_masks.I, label="Маски (β↓50%): R0=$(round(R0_masks, digits=2))")
plot!(plt_control, df_treatment.t, df_treatment.I, label="Лечение (γ↑50%): R0=$(round(R0_treatment, digits=2))")
plt_control

# **Вывод:** Все меры эффективны, но карантин (снижение контактов) 
# дает самый быстрый эффект в краткосрочной перспективе.

# # Сохранение графиков

plots_path = joinpath(@__DIR__, "..", "plots", "parameter_study_literate")
mkpath(plots_path)

savefig(plt_β, joinpath(plots_path, "beta_study.png"))
savefig(plt_c, joinpath(plots_path, "c_study.png"))
savefig(plt_γ, joinpath(plots_path, "gamma_study.png"))
savefig(plt_control, joinpath(plots_path, "control_measures.png"))

println("\nВсе графики сохранены в: ", plots_path)

# # Сводная таблица результатов

println("\n" * "="^60)
println("СВОДНАЯ ТАБЛИЦА РЕЗУЛЬТАТОВ")
println("="^60)

results = []
for β in β_values
    p = [β, c_fixed, γ_fixed]
    df, R0 = run_sir(p)
    peak_idx = argmax(df.I)
    push!(results, (Параметр="β=$(β)", R0=round(R0, digits=2), 
                   Пик=round(df.I[peak_idx], digits=1), 
                   Время_пика=round(df.t[peak_idx], digits=1),
                   Итоговые_R=round(df.R[end], digits=1)))
end

for c in c_values
    p = [β_fixed, c, γ_fixed]
    df, R0 = run_sir(p)
    peak_idx = argmax(df.I)
    push!(results, (Параметр="c=$(c)", R0=round(R0, digits=2), 
                   Пик=round(df.I[peak_idx], digits=1), 
                   Время_пика=round(df.t[peak_idx], digits=1),
                   Итоговые_R=round(df.R[end], digits=1)))
end

for γ in γ_values
    p = [β_fixed, c_fixed, γ]
    df, R0 = run_sir(p)
    peak_idx = argmax(df.I)
    push!(results, (Параметр="γ=$(γ)", R0=round(R0, digits=2), 
                   Пик=round(df.I[peak_idx], digits=1), 
                   Время_пика=round(df.t[peak_idx], digits=1),
                   Итоговые_R=round(df.R[end], digits=1)))
end

# Выводим таблицу
for r in results
    println(r.Параметр, ": R0=", r.R0, ", пик=", r.Пик, " на ", r.Время_пика, " день, итого переболело=", r.Итоговые_R)
end

# ## Заключение
# 
# Проведенное исследование показывает, что:
# 1. **R0** является ключевым параметром, определяющим динамику эпидемии
# 2. Все три параметра (β, c, γ) линейно влияют на R0
# 3. Для контроля эпидемии наиболее эффективно снижение числа контактов (карантин)
#   и использование защитных мер (маски)
# 4. Комбинация мер дает наилучший результат

println("\nИсследование параметров модели SIR завершено!")

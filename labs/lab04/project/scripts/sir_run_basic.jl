# ## Базовый эксперимент SIR-модели
# 
# ### Описание
# 
# Данный скрипт запускает базовый эксперимент с фиксированными параметрами 
# для проверки работоспособности агентной SIR-модели.

using DrWatson
@quickactivate "project"
using Agents, DataFrames, Plots, JLD2

include(srcdir("sir_model.jl"))

# Параметры эксперимента
params = Dict(
    :Ns => [1000, 1000, 1000],
    :β_und => [0.5, 0.5, 0.5],
    :β_det => [0.05, 0.05, 0.05],
    :infection_period => 14,
    :detection_time => 7,
    :death_rate => 0.02,
    :reinfection_probability => 0.1,
    :Is => [0, 0, 1],
    :seed => 42,
    :n_steps => 100,
)

# Инициализация модели
# Для initialize_sir берем только нужные параметры
model_params = Dict(k => v for (k, v) in params if k != :n_steps)
model = initialize_sir(; model_params...)

# Подготовка массивов для хранения данных
times = Int[]
S_vals = Int[]
I_vals = Int[]
R_vals = Int[]
total_vals = Int[]

# Запуск симуляции вручную
for step in 1:params[:n_steps]
    Agents.step!(model, 1)
    
    push!(times, step)
    push!(S_vals, susceptible_count(model))
    push!(I_vals, infected_count(model))
    push!(R_vals, recovered_count(model))
    push!(total_vals, total_count(model))
end

# Создаём DataFrame
agent_df = DataFrame(time = times, susceptible = S_vals, infected = I_vals, recovered = R_vals)
model_df = DataFrame(time = times, total = total_vals)

# Визуализация
plot(agent_df.time, agent_df.susceptible, label = "Восприимчивые", xlabel = "Дни", ylabel = "Количество")
plot!(agent_df.time, agent_df.infected, label = "Инфицированные")
plot!(agent_df.time, agent_df.recovered, label = "Выздоровевшие")
plot!(agent_df.time, model_df.total, label = "Всего (включая умерших)", linestyle = :dash)
savefig(plotsdir("sir_basic_dynamics.png"))

# Сохранение данных
@save datadir("sir_basic_agent.jld2") agent_df
@save datadir("sir_basic_model.jld2") model_df

# ### Выводы

# 1. Модель успешно инициализирована и запущена. Динамика численности агентов
#    соответствует классической SIR-модели: восприимчивые снижаются,
#    инфицированные проходят через пик, выздоровевшие нарастают.
#
# 2. При заданных параметрах ($\beta > \gamma$) наблюдается развитие эпидемии.
#    Репродуктивное число $R_0 = \beta / \gamma$ превышает пороговое значение 1,
#    что подтверждается динамикой модели.
#
# 3. В ходе эпидемического процесса формируется популяционный иммунитет:
#    значительная часть восприимчивых агентов переходит в группу выздоровевших.
#
# 4. Наличие летальных исходов приводит к снижению общей численности популяции,
#    что отражается на графике пунктирной линией.
#
# 5. Полученные данные и график сохранены для дальнейшего анализа и сравнения
#    с другими сценариями (сканирование параметров, влияние миграции и т.д.).

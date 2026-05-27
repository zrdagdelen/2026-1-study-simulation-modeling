#!/usr/bin/env julia
"""
sir_model.jl - Ядро дискретно-событийной SIR-модели
Автор: Лабораторная работа №8
Дата: 2026
"""

using ResumableFunctions
using ConcurrentSim
using Distributions
using DataFrames
using Random

# ============================================================
# Вспомогательные функции для обновления массивов статистики
# ============================================================

"""
    increment!(a::Array{Int64})

Добавляет в конец массива значение, увеличенное на 1 относительно последнего элемента.
Используется для фиксации увеличения численности группы.
"""
function increment!(a::Array{Int64})
    push!(a, a[end] + 1)
end

"""
    decrement!(a::Array{Int64})

Добавляет в конец массива значение, уменьшенное на 1 относительно последнего элемента.
Используется для фиксации уменьшения численности группы.
"""
function decrement!(a::Array{Int64})
    push!(a, a[end] - 1)
end

"""
    carryover!(a::Array{Int64})

Дублирует последнее значение массива.
Используется для фиксации отсутствия изменений в численности группы.
"""
function carryover!(a::Array{Int64})
    push!(a, a[end])
end

# ============================================================
# Определение структур данных
# ============================================================

"""
    mutable struct SIRPerson

Агент-индивид в популяции.
- id::Int64: уникальный идентификатор
- status::Symbol: текущее состояние (:S, :I, :R)
"""
mutable struct SIRPerson
    id::Int64
    status::Symbol  # :S, :I, :R
end

"""
    mutable struct SIRModel

Хранит всё состояние модели.
- sim::ConcurrentSim.Simulation: объект симуляции
- β::Float64: вероятность заражения при контакте
- c::Float64: частота контактов
- γ::Float64: интенсивность выздоровления
- μ::Float64: интенсивность смертности (опционально)
- ta::Array{Float64}: временные метки событий
- Sa::Array{Int64}: временной ряд восприимчивых
- Ia::Array{Int64}: временной ряд инфицированных
- Ra::Array{Int64}: временной ряд переболевших
- allIndividuals::Array{SIRPerson}: массив всех индивидов
"""
mutable struct SIRModel
    sim::ConcurrentSim.Simulation
    β::Float64
    c::Float64
    γ::Float64
    μ::Float64  # дополнительный параметр для демографии
    ta::Array{Float64}
    Sa::Array{Int64}
    Ia::Array{Int64}
    Ra::Array{Int64}
    allIndividuals::Array{SIRPerson}
end

# ============================================================
# Функции обновления статистики при событиях
# ============================================================

"""
    infection_update!(sim::ConcurrentSim.Simulation, m::SIRModel)

Вызывается в момент заражения. Фиксирует текущее время,
уменьшает число восприимчивых, увеличивает число инфицированных.
"""
function infection_update!(sim::ConcurrentSim.Simulation, m::SIRModel)
    push!(m.ta, ConcurrentSim.now(sim))
    decrement!(m.Sa)
    increment!(m.Ia)
    carryover!(m.Ra)
end

"""
    recovery_update!(sim::ConcurrentSim.Simulation, m::SIRModel)

Вызывается при выздоровлении. Фиксирует текущее время,
уменьшает число инфицированных, увеличивает число переболевших.
"""
function recovery_update!(sim::ConcurrentSim.Simulation, m::SIRModel)
    push!(m.ta, ConcurrentSim.now(sim))
    carryover!(m.Sa)
    decrement!(m.Ia)
    increment!(m.Ra)
end

# ============================================================
# Основной процесс агента
# ============================================================

"""
    live(env::ConcurrentSim.Simulation, individual::SIRPerson, m::SIRModel)

Основная логика жизненного цикла индивида.
Использует @resumable для приостановки выполнения без блокировки.
"""
@resumable function live(env::ConcurrentSim.Simulation, 
                         individual::SIRPerson, 
                         m::SIRModel)
    while individual.status == :S
        # Ожидание до следующего контакта
        @yield timeout(env, rand(Exponential(1/m.c)))
        
        # Выбор случайного собеседника
        alter = individual
        while alter == individual
            N = length(m.allIndividuals)
            index = rand(DiscreteUniform(1, N))
            alter = m.allIndividuals[index]
        end
        
        # Попытка заражения
        if alter.status == :I
            if rand(Uniform(0, 1)) < m.β
                individual.status = :I
                infection_update!(env, m)
            end
        end
    end
    
    # Цикл для инфицированных
    if individual.status == :I
        # Ожидание выздоровления (стохастическое)
        @yield timeout(env, rand(Exponential(1/m.γ)))
        individual.status = :R
        recovery_update!(env, m)
    end
end

# ============================================================
# Функции управления моделью
# ============================================================

"""
    MakeSIRModel(uθ, p; μ=0.0)

Создаёт экземпляр SIRModel.
- uθ: начальные условия [S₀, I₀, R₀]
- p: параметры [β, c, γ]
- μ: интенсивность смертности (по умолчанию 0)
"""
function MakeSIRModel(uθ, p; μ=0.0)
    (S, I, R) = uθ
    N = S + I + R
    (β, c, γ) = p
    
    sim = ConcurrentSim.Simulation()
    
    # Создание массива индивидов
    allIndividuals = SIRPerson[]
    for i in 1:S
        push!(allIndividuals, SIRPerson(i, :S))
    end
    for i in (S+1):(S+I)
        push!(allIndividuals, SIRPerson(i, :I))
    end
    for i in (S+I+1):N
        push!(allIndividuals, SIRPerson(i, :R))
    end
    
    # Инициализация временных рядов
    ta = Float64[0.0]
    Sa = Int64[S]
    Ia = Int64[I]
    Ra = Int64[R]
    
    SIRModel(sim, β, c, γ, μ, ta, Sa, Ia, Ra, allIndividuals)
end

"""
    activate(m::SIRModel)

Запускает все процессы агентов. Процессы регистрируются в симуляции.
"""
function activate(m::SIRModel)
    [@process live(m.sim, individual, m) for individual in m.allIndividuals]
end

"""
    sir_run(m::SIRModel, tf::Float64)

Запускает симуляцию до конечного времени tf.
"""
function sir_run(m::SIRModel, tf::Float64)
    ConcurrentSim.run(m.sim, tf)
end

"""
    out(m::SIRModel) -> DataFrame

Собирает результаты в DataFrame с колонками :t, :S, :I, :R.
"""
function out(m::SIRModel)
    result = DataFrame()
    result[:, :t] = m.ta
    result[:, :S] = m.Sa
    result[:, :I] = m.Ia
    result[:, :R] = m.Ra
    return result
end

# ============================================================
# Функции для дополнительных заданий
# ============================================================

"""
    MakeSIRModel_deterministic_recovery(uθ, p)

Создаёт модель с детерминированным временем выздоровления.
"""
function MakeSIRModel_deterministic_recovery(uθ, p; μ=0.0)
    (S, I, R) = uθ
    N = S + I + R
    (β, c, γ) = p
    
    sim = ConcurrentSim.Simulation()
    
    allIndividuals = SIRPerson[]
    for i in 1:S
        push!(allIndividuals, SIRPerson(i, :S))
    end
    for i in (S+1):(S+I)
        push!(allIndividuals, SIRPerson(i, :I))
    end
    for i in (S+I+1):N
        push!(allIndividuals, SIRPerson(i, :R))
    end
    
    ta = Float64[0.0]
    Sa = Int64[S]
    Ia = Int64[I]
    Ra = Int64[R]
    
    SIRModel(sim, β, c, γ, μ, ta, Sa, Ia, Ra, allIndividuals)
end

"""
    live_deterministic(env, individual, m)

Версия live с детерминированным временем выздоровления.
"""
@resumable function live_deterministic(env::ConcurrentSim.Simulation,
                                        individual::SIRPerson,
                                        m::SIRModel)
    while individual.status == :S
        @yield timeout(env, rand(Exponential(1/m.c)))
        
        alter = individual
        while alter == individual
            N = length(m.allIndividuals)
            index = rand(DiscreteUniform(1, N))
            alter = m.allIndividuals[index]
        end
        
        if alter.status == :I
            if rand(Uniform(0, 1)) < m.β
                individual.status = :I
                infection_update!(env, m)
            end
        end
    end
    
    if individual.status == :I
        # Детерминированное время выздоровления
        @yield timeout(env, 1/m.γ)
        individual.status = :R
        recovery_update!(env, m)
    end
end

"""
    activate_deterministic(m::SIRModel)

Активация с детерминированным выздоровлением.
"""
function activate_deterministic(m::SIRModel)
    [@process live_deterministic(m.sim, individual, m) for individual in m.allIndividuals]
end

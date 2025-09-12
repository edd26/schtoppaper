"DataSource.jl" |> srcdir |> include
"""
Abstract type to bind signal based sources.
"""
abstract type SignalSource <: DataSource end

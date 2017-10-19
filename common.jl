# Commomly used functions

#FIXME: Comparing fp numbers using fixed epsilion is prone to error
equal(a::dp, b::dp, epsilon::dp = PRECIS) =  abs(a - b) < epsilon
locate(item, list::Vector) = findfirst(x -> x == item, list)
locate(fp::dp, list::Vector{dp}, epsilon::dp = PRECIS) = findfirst(x -> equal(x, fp, epsilon), list)

# Atomic symbols and numbers
atomic_sym(at_num::Int) = p_table[at_num]
atomic_num(at_sym::T) where {T<:AbstractString} = locate(at_sym, p_table)
atomic_sym(at_num::Vector{Int}) = [atomic_sym(at_num[i]) for i in 1:length(at_num)]
atomic_num(at_sym::Vector{T}) where {T<:AbstractString} = [atomic_num(at_sym[i]) for i in 1:length(at_sym)]

# Helper functions for filenames
wd_fname(wd::AbstractString, fname::AbstractString) = endswith(wd, "/") ? wd * fname : wd * "/" * fname
wd_tagfname(wd::AbstractString, tag::Int, prefix::AbstractString = "") = endswith(wd, "/") ? wd * prefix * "$tag" : wd * "/" * prefix * "$tag"

# Parsing
parse_list(dt::DataType, str::Union{Vector{T}, SubArray{T, 1}}) where {T<:AbstractString} = [parse(dt, str[i]) for i in 1:length(str)]

# Distance between point a and b
dist(a::Vector{dp}, b::Vector{dp}) = norm(a - b)

# Unit cell related
# Reciprocal lattice in 2π unit
reciprocal(lattice::Matrix{dp}) = inv(lattice).'
# Volume
volume(lattice::Matrix{dp}) = det(lattice)


# Simple polynomial fitting
# increasing order
polyfit_inc(x::Vector{dp}, y::Vector{dp}, order::Int) = [x[i]^p for i = 1:length(x), p = 0:order] \ y
# decreasing order
polyfit_dec(x::Vector{dp}, y::Vector{dp}, order::Int) = M = [x[i]^p for i = 1:length(x), p = order:-1:0] \ y
# default is lowest orders first (increasing order)
polyfit(x::Vector{dp}, y::Vector{dp}, order::Int) = polyfit_inc(x, y, order)

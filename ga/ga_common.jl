prob(p::Real) = rand() < p

swap!(v::Vector, i::Int, j::Int) = v[i], v[j] = v[j], v[i]

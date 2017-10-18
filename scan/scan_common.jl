# from Julia examples/ndgrid.jl
function meshgrid(vx::AbstractVector{T}, vy::AbstractVector{T}) where T
    m, n = length(vy), length(vx)
    vx = reshape(vx, 1, n)
    vy = reshape(vy, m, 1)
    (repmat(vx, m, 1), repmat(vy, 1, n))
end

function ndgrid(v1::AbstractVector{T}, v2::AbstractVector{T}) where T
    m, n = length(v1), length(v2)
    v1 = reshape(v1, m, 1)
    v2 = reshape(v2, 1, n)
    (repmat(v1, 1, n), repmat(v2, m, 1))
end

function move!(geo::Geometry, at::AtomList, vect::Vector{dp})

    for i = 1:length(at)
        for j = 1:3
            geo.coord[j, at[i]] += vect[j]
        end
    end

    return geo
end

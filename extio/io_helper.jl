# Helper functions for filenames
wd_fname(wd::AbstractString, fname::AbstractString) = endswith(wd, "/") ? wd * fname : wd * "/" * fname
wd_tagfname(wd::AbstractString, tag::Int, prefix::AbstractString = "") = endswith(wd, "/") ? wd * prefix * "$tag" : wd * "/" * prefix * "$tag"

# Batch copying files
function batch_cp(flist::Vector{T}, src::T, dst::T) where {T<:AbstractString}

    for i = 1:length(flist)
        fsrc = wd_fname(src, flist[i])
        fdst = wd_fname(dst, flist[i])
        cp(fsrc, fdst, remove_destination = true)
    end

end

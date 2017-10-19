mutable struct ScanWd

    wd::String
    subwd::Vector{ScanWd}

    function ScanWd(wd::AbstractString, rsc::RepScale)

        swd = new()
        swd.wd = string(wd)
        swd.subwd = ScanWd[]

        for i = 1:length(rsc.scaled)
            new_wd = wd_tagfname(wd, i, SCALE_PREFIX)
            new_swd = ScanWd(new_wd, rsc.scaled[i])
            push!(swd.subwd, new_swd)
        end

        return swd
    end

end

function flatten(swd::ScanWd)

    f = String[]
    push!(f, swd.wd)
    for i = 1:length(swd.subwd)
        f = vcat(f, flatten(swd.subwd[i]))
    end

    return f
end

module SWSTimes

using Base.Dates

using SWSs

export
    SWSTime,
    read_sheba_time

immutable SWSTime <: AbstractSWS
    s::SWS
    time::DateTime
end

const SWST = Union{SWSTime, Array{SWSTime}}

"""
    getindex(s::Union{SWSTime,Array{SWSTime}}, k::Symbol) -> v

Return an array of values `v` for the fields named `k` in an array of `SWSTime`s.

Use `:time` for `k` to get the `DateTime`s for each measurement; otherwise, use `:sws`
to get an array of `SWS`s; or finally use any fieldname in the `SWS` type to get that.
"""
Base.getindex(s::Array{SWSTime}, k::Symbol) = k == :time ? [ss.time for ss in s] :
                                                  k == :sws ? [ss.s for ss in s] :
                                                  [getfield(ss.s, k) for ss in s]
Base.getindex(s::SWSTime, k::Symbol) = k == :time ? s.time : s.s[k]

"""
    read_sheba_time(file, args...; comment=r\"^[#%]\") -> st::Vector{SWSTime}

Read results from a SHEBA-formatted file, including the date and time of the
events listed therein, returning the values `st`.

Optionally specify the Regex used to determine which lines are comments with `comment`.
"""
function read_sheba_time(file, args...; comment=r"^[#%]")
    s = read_sheba(file, args...)
    n = length(s)
    st = Vector{SWSTime}(n)
    data = readdlm(file)
    data = data[.!ismatch.(comment, string.(data[:,1])), :]
    size(data, 1) == n || error("Length of splits read ($n) does "*
        "not match the number of non-comment lines ($(size(data, 1))) in file '$file'")
    for i in 1:n
        y = data[i,1]÷1000
        d = data[i,1]%1000
        h = data[i,2]÷100
        m = data[i,2]%100
        date = Date(y) + Day(d)
        st[i] = SWSTime(s[i], DateTime(year(date), month(date), day(date), h, m))
    end
    st
end

end # module

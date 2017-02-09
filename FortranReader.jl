"""
# Module FortranReader

FortranReader provides basic commands to read the 'unformatted' files written
by Fortran programs.

Although not technically portable, most Fortran programs write these files in
a predictable way, in 'records'.  Each record contains either one or several
variables, marked before and after with the length of the variable(s) in bytes.
This marker is written as a 4-byte integer.

function `read_record(::IO, ::DataType, rank::Int, dims)` returns the value in
the next record from the stream and data type given.
"""
module FortranReader

export
    read_record,
    skip_record

# Standard Fortran number sizes
typealias Real4 Float32
typealias Real8 Float64
typealias Integer4 Int32
typealias Integer8 Int64
typealias Complex4 Complex64
typealias Complex8 Complex128
typealias Character String

const len4 = sizeof(Integer4)
const len8 = sizeof(Integer8)

"""
    read_record(f::IO, T::DataType, dims...) -> v::T

Return the next record in the stream `f`, in the format of `T`.  This means,
`read_record` assumes that the type of `v` is `T` (e.g., `Float64`).
`dims...` holds the dimensions; the default (no argument) assumes a scalar record.
"""
function read_record(f::IOStream, T::DataType, dims...)
    all(Bool[typeof(dim) <: Integer for dim in dims]) || error("`dims` must be integers")
    d, len = read_raw(f)
    Tlen = sizeof(T)
    if T <: String
        return String(d)
    end
    rank = length(dims)
    if rank == 0
        len == Tlen || error("Record at $(position(f) - len4 - len) is not correct " *
            "length for scalar of type $T")
        reinterpret(T, d)[1]
    elseif rank == 1
        n = len ÷ Tlen
        n != dims[1] && warn("Length of record does not match dimension supplied " *
            "(Requested $(dims[1]); read size $n)")
        reinterpret(T, d)
    elseif rank >= 2
        len ÷ Tlen == prod(dims) || error("Requested dimensions do not match size "*
            "of record.  (Requested $dims (size $(prod(dims))); record size $(len÷Tlen))")
        reshape(reinterpret(T, d), dims)
    else
        error("Arrays must have positive rank")
    end
end

"""
    skip_record(f::IOStream) -> n::Int

Skip over the next record, returning the length of the record.
"""
function skip_record(f::IOStream)
    nstart = reinterpret(Integer4, read(f, len4))[1]
    pos = position(f)
    skip(f, nstart)
    nend = reinterpret(Integer4, read(f, len4))[1]
    nstart == nend || error("Error skipping over record at byte $pos " *
        "in stream $(f.name)")
    nstart
end

"""
    read_raw(f) -> d::Array{UInt8,1}, n::Int

Return the raw data from stream `f`, checked that the start and end records match,
and the number of bytes `n` read.
"""
function read_raw(f::IOStream)
    nstart = reinterpret(Integer4, read(f, len4))[1]
    d = read(f, nstart)
    nend = reinterpret(Integer4, read(f, len4))[1]
    nstart == nend || error("Error reading value of type $T from stream $(f.name)")
    d, nstart
end

end # module

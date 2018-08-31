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
    skip_record,
    write_record

# Standard Fortran number sizes
const Real4 = Float32
const Real8 = Float64
const Integer4 = Int32
const Integer8 = Int64
const Complex4 = Complex{Float32}
const Complex8 = Complex{Float64}
const Character = String

const len4 = sizeof(Integer4)
const len8 = sizeof(Integer8)

const permissible_types = (Real4, Real8, Integer4, Integer8, Complex4, Complex8, Character)

"""
    read_record(f::IO, T::DataType, dims...) -> v::T

Return the next record in the stream `f`, in the format of `T`.  This means,
`read_record` assumes that the type of `v` is `T` (e.g., `Float64`).
`dims...` holds the dimensions; the default (no argument) assumes a scalar record.
"""
function read_record(f::IOStream, T::DataType, dims...)
    all(Bool[typeof(dim) <: Integer for dim in dims]) || error("`dims` must be integers")
    d, len = read_raw(f)
    if T <: String
        return String(d)
    end
    Tlen = sizeof(T)
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
    pos = position(f)
    d = read(f, nstart)
    nend = reinterpret(Integer4, read(f, len4))[1]
    nstart == nend || error("Error reading raw data at position $pos from stream "
        * "$(f.name): record length markers do not match ($nstart != $nend)")
    d, nstart
end

"""
    write_record(f::IOStream, T::DataType, x...) -> n::Int

Write a record to the stream `f` containing `x`, which will be converted to have
element type `T`.  `T` must therefore correspond to one of the Fortran types 
`Float{32,64}`, `Int{32,64}`, `Complex{64,128}` or `String`.  You can also use the
equivalent Fortran names which are exported by the module:  `Integer{4,8}`, `Real{4,8}`,
`Complex{4,8}` and `Character`.

Return the **total** number of bytes written to the stream `f`, including the start-
and end-markers of the record.
"""
function write_record(f::IOStream, T::DataType, xs...)
    T in permissible_types || error("Type of record $T is not supported.  " *
        "Choose one of $permissible_types")
    bytes_written = 0
    for x in xs
        n = T == Character ? Integer4(length(x)) : Integer4(sizeof(eltype(T))*length(x))
        write(f, n)
        if T == Character
            write(f, x)
        elseif ndims(x) > 0
            write(f, convert(Array{T}, x))
        else
            write(f, convert(T, x))
        end
        write(f, n)
        bytes_written += 2*len4 + n
    end
    bytes_written
end

end # module

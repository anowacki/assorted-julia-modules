"""
# TarFile

TarFile provides methods to extract files in a tar archive into memory.

The main function is `extract`, which returns an array of UInt8 arrays,
each corresponding to each file in the archive, and a file name for each.

## Dependencies
This module relies on the python module `tarfile` to do all the work and simply
provides a convenient wrapper.  Hence you need:

- PyCall
- A python installation with `tarfile`.

## Using
```
julia> using TarFile

julia> content, names = extract("/path/to/tarfile.tar")
```

## Limitations
- TarFile only operates on tar files on disk at present.
"""
module TarFile

using PyCall

export extract

@pyimport tarfile as TF

"""
    extract(file) -> file_contents::Array{Array{Uint8}}, file_names::Array{String}

Read a tape archive ('tar') file and return the contents as `file_contents` and
file paths as `file_names` to all the members of the archive.

Archives with gzip or bzip2 compression are supported without any further action
required by the user.
"""
function extract(file)
    # If this gets too slow, look at:
    #   https://blogs.it.ox.ac.uk/inapickle/2011/06/20/high-memory-usage-when-using-pythons-tarfile-module/
    f = TF.open(file)
    contents = Array{Array{UInt8,1},1}()
    names = Array{String,1}()
    for name in f[:getnames]()
        push!(names, name)
        push!(contents, f[:extractfile](name)[:read]())
    end
    f[:close]()
    contents, names
end

end # module

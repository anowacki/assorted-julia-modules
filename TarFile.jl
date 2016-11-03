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

import Glob: FilenameMatch, @fn_str
using PyCall

export extract

@pyimport tarfile as TF

"""
    extract(file, pattern="*", verbose=false) -> file_contents::Array{Array{Uint8}}, file_names::Array{String}

Read a tape archive ('tar') file and return the contents as `file_contents` and
file paths as `file_names` to all the members of the archive.  If `pattern` is
specified, then only return files matching the pattern.  This can be a Regex
(specified by `r"..."`) or a file globbing pattern.

Archives with gzip or bzip2 compression are supported without any further action
required by the user.
"""
function extract(file, pattern::Union{Regex,FilenameMatch}=r""; verbose=false)
    # If this gets too slow, look at:
    #   https://blogs.it.ox.ac.uk/inapickle/2011/06/20/high-memory-usage-when-using-pythons-tarfile-module/
    verbose && info("Opening tar file '$file'")
    contents = Array{Array{UInt8,1},1}()
    names = Array{String,1}()
    f = TF.open(file)
    try
        for member in f
            filename = member[:path]
            if ismatch(pattern, filename) && member[:isfile]()
                handle = f[:extractfile](member)
                handle == nothing && continue
                verbose && info("  Extracting file $filename")
                push!(names, filename)
                push!(contents, handle[:read]())
            end
        end
    catch err
        verbose && info("  Cleaning up after encountering error:\n$err")
    finally
        f[:close]()
    end
    contents, names
end
extract(file, pattern::AbstractString; kw...) = extract(file, FilenameMatch(pattern); kw...)

end # module

# B+ and DataStructures both define `update!`
const update! = Bplus.update!

"
Retrieves the key from an ordered dictionary (or element from an ordered set),
    by its index.

The OrderedCollections package doesnt have an official way to do this yet;
    see https://github.com/JuliaCollections/OrderedCollections.jl/issues/93.
"
function ordered_collection_get_idx(oc::Union{OrderedDict, OrderedSet}, idx::Integer)
    # Make sure cached ops have all been applied.
    iterate(oc)

    inner = if oc isa OrderedDict
        oc
    elseif oc isa OrderedSet
        oc.dict
    else
        error("Unhandled: ", typeof(oc))
    end

    return inner.keys[idx]
end

"Scoped wrapper for read-write locks"
@inline function lock_read(to_do, rwl::ReadWriteLock)
    readlock(rwl)
    try
        return to_do()
    finally
        readunlock(rwl)
    end
end
"Scoped wrapper for read-write locks"
@inline function lock_write(to_do, rwl::ReadWriteLock)
    lock(rwl)
    try
        return to_do()
    finally
        unlock(rwl)
    end
end
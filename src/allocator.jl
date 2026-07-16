
"
Some kind of memory allocator, mainly used for the grids.
Default behavior of its interface is to allocate things from the heap
   and let the GC handle them on release.
"
abstract type AbstractMarkovAllocator end
"Pass an empty tuple for the size in order to get an empty 1D vector"
function markov_allocator_acquire_array(allocator::AbstractMarkovAllocator,
                                        size::Tuple{Vararg{Integer}},
                                        ::Type{T}
                                       )::AbstractArray{T, (isempty(size) ? 1 : length(size))} where {T}
    if isempty(size)
        return preallocated_vector(T, 128)
    else
        return Array{T, length(size)}(undef, size)
    end
end
markov_allocator_acquire_ordered_set(allocator::AbstractMarkovAllocator,
                                     T::Type
                                    )::OrderedSet{T} = OrderedSet{T}()
markov_allocator_acquire_set(allocator::AbstractMarkovAllocator, T::Type)::Set{T} = Set{T}()
markov_allocator_release_array(allocator::AbstractMarkovAllocator, data::AbstractArray) = nothing
markov_allocator_release_ordered_set(allocator::AbstractMarkovAllocator, s::OrderedSet) = nothing
markov_allocator_release_set(allocator::AbstractMarkovAllocator, s::Set) = nothing
Base.close(::AbstractMarkovAllocator) = nothing


########

"Allocates on the heap and lets Julia's GC take the array afterwards"
struct MarkovAllocatorHeap <: AbstractMarkovAllocator end


"Allocates on the heap and re-uses instances as much as possible"
struct MarkovAllocatorHeapReused <: AbstractMarkovAllocator
    arrays::Dict{Tuple{DataType, Int}, Vector{<:Array}}
    ordered_sets::Dict{DataType, Vector{<:OrderedSet}}
    sets::Dict{DataType, Vector{<:Set}}
end
MarkovAllocatorHeapReused() = MarkovAllocatorHeapReused(
    Dict{Tuple{DataType, Int}, Vector{<:Array}}(),
    Dict{DataType, Vector{<:OrderedSet}}(),
    Dict{DataType, Vector{<:Set}}()
)

"Helper struct that exists to return a specific value for 'Base.length()'"
struct ConstantLength
    n::Int
end
Base.length(c::ConstantLength) = c.n

WARN_REUSED_HEAP_ARRAY_IMPL_CHANGED::Bool = false

function markov_allocator_acquire_array(alloc::MarkovAllocatorHeapReused,
                                        desired_size::TSize,
                                        ::Type{TElement}
                                       )::Array{TElement} where {TSize<:Tuple{Vararg{Integer}}, TElement}
    N = isempty(desired_size) ? 1 : length(desired_size)

    existing::Vector{Array{TElement, N}} = get!(
        () -> preallocated_vector(Array{TElement, N}, 128),
        alloc.arrays, (TElement, N)
    )

    # Find the first array that matches the size.
    # For faster searching, arrays are always sorted by their total length.
    if isempty(desired_size)
        if isempty(existing)
            return preallocated_vector(TElement, 128)
        else
            # For now I think the best solution is to always return the last/largest array.
            # I think the cost of deleting from the center/back of a vector
            #    outweights the impact of giving a large presized array to a small user.
            v = pop!(existing)
            empty!(v)
            return v
        end
    else
        n_elements = prod(desired_size)

        i = searchsortedfirst(existing, ConstantLength(n_elements), by=length)
        while (i <= length(existing)) && (length(existing[i]) == n_elements) && (size(existing[i]) != desired_size)
            i += 1
        end

        if (i <= length(existing)) && (size(existing[i]) == desired_size)
            a = existing[i]
            deleteat!(existing, i)
            return a
        else
            # We could try to resize a larger vector, but the user wanted a specific size
            #    so that extra memory is more likely to be wasted here and needed elsewhere.
            return Array{TElement, N}(undef, desired_size)
        end
    end
end
function markov_allocator_release_array(alloc::MarkovAllocatorHeapReused,
                                        array::Array{TElement, N}) where {TElement, N}
    # For vectors, resize them to max out their underlying memory before storing them.
    # This unfortunately requires some internal details of Array.
    if (N == 1)
        if hasfield(typeof(array), :ref) && hasfield(typeof(array.ref), :mem)
            capacity = length(array.ref.mem)
            resize!(array, capacity)
        elseif !WARN_REUSED_HEAP_ARRAY_IMPL_CHANGED
            WARN_REUSED_HEAP_ARRAY_IMPL_CHANGED = true
            @warn "Internal implementation details of `Array` have changed, and we lost some allocator efficiency! Find this message for more details."
        end
    end

    # Keep the array lists sorted by total length.
    existing::Vector{Array{TElement, N}} = alloc.arrays[(TElement, N)]
    insert!(
        existing,
        1 + searchsortedlast(existing, ConstantLength(length(array)), by=length),
        array
    )

    return nothing
end

function markov_allocator_acquire_ordered_set(alloc::MarkovAllocatorHeapReused, ::Type{T}) where {T}
    existing::Vector{OrderedSet{T}} = get!(
        () -> preallocated_vector(OrderedSet{T}, 16),
        alloc.ordered_sets, T
    )
    if isempty(existing)
        return OrderedSet{T}()
    else
        @markovjunior_assert(isempty(existing[end]))
        return pop!(existing)
    end
end
function markov_allocator_release_ordered_set(alloc::MarkovAllocatorHeapReused, os::OrderedSet{T}) where {T}
    empty!(os)
    push!(alloc.ordered_sets[T]::Vector{OrderedSet{T}}, os)
    return nothing
end

function markov_allocator_acquire_set(alloc::MarkovAllocatorHeapReused, ::Type{T}) where {T}
    existing::Vector{Set{T}} = get!(
        () -> preallocated_vector(Set{T}, 16),
        alloc.sets, T
    )
    if isempty(existing)
        return Set{T}()
    else
        @markovjunior_assert(isempty(existing[end]))
        return pop!(existing)
    end
end
function markov_allocator_release_set(alloc::MarkovAllocatorHeapReused, s::Set{T}) where {T}
    empty!(s)
    push!(alloc.sets[T]::Vector{Set{T}}, s)
    return nothing
end
function weighted_random_array_element(weights_iter, total_sum_weight, uniform_random_float)::Int
    F = typeof(uniform_random_float)
    if !isempty(weights_iter)
        F = promote_type(F, typeof(first(weights_iter)))
    end

    threshold = convert(F, uniform_random_float * total_sum_weight)
    current_level = zero(F)
    max_i = 0
    for (i, weight) in enumerate(weights_iter)
        max_i = i
        if current_level + weight >= threshold
            return i
        else
            current_level += weight
        end
    end

    # If we get here, then our random float was near 1.0
    #    and floating-point error kept us from choosing the last element.
    return max_i
end


"
A purely-cosmetic type annotation for readability; actually turns into `Any`.
Example usage:

````
columns::@IterOf{AbstractVector{Int}} = (@view m[:, i] for i in 2:6)
````
"
macro IterOf(expr_brace)
    @bp_check(Base.isexpr(expr_brace, :braces))
    return Any
end


@inline function is_symmetric_along(array::AbstractArray, _axis::Integer,   element_comparator = isequal)::Bool
    axis = convert(Int, _axis)
    if axis > ndims(array)
        # In higher dimensions this array only has one slice, so it's symmetrical!
        return true
    end

    N = ndims(array)
    sliced(i) = view(array, ntuple(j-> (j==Axis) ? i : Colon(), Val(N))...)

    n = size(array, axis)
    for i in 1:(n÷2)
        i2 = n - i + 1
        element_pairs = zip(selectdim(array, axis, i), selectdim(array, axis, i2))
        if any(!element_comparator(a, b) for (a,b) in element_pairs)
            return false
        end
    end

    return true
end

@inline function are_axes_swappable(array::AbstractArray, axes::NTuple{2, Integer}, flip::Bool, element_comparator = isequal)::Bool
    A1, A2 = convert.(Ref(Int), axes)
    N = ndims(array)
    if (A1 > N) || (A2 > N)
        return (A1 > N) && (A2 > N)
    elseif size(array, A1) != size(array, A2)
        return false
    elseif size(array, A1) == 0
        return true
    end

    if size(selectdim(array, A1, 1)) != size(selectdim(array, A2, 1))
        return false
    end
    compare_slices(i, j) = let (s1,s2) = (selectdim(array, A1, i), selectdim(array, A2, j))
        @markovjunior_assert(size(s1) == size(s2),
                             N, "/", A1, "-", A2, "/", i, "-", j, "/", size(s1), "-", size(s2))
        all(element_comparator(a, b) for (a, b) in zip(s1, s2))
    end

    n = size(array, A1)
    return if flip
        all(compare_slices(i, n - i + 1) for i in 1:n)
    else
        all(compare_slices(i, i) for i in 1:n)
    end
end
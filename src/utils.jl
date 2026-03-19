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
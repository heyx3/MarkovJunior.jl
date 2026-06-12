# The C API for this package.
# Internal stuff is prefixed "lib_"
# External stuff is prefixed "jmj_"

const Lib_ID = Cuint

# Parsed algorithms are stored in a global lookup by unique ID.
# Otherwise they'd be a raw pointer and the GC could take them in between library calls!
const LIB_ALGORITHMS = Dict{Lib_ID, MarkovAlgorithm}()
const LIB_NEXT_ALGORITHM_ID = Ref{Lib_ID}(1)
# There's no telling how many threads will access this package,
#    so we need locks around the global store.
const LIB_ALGORITHM_LOCKER = ReadWriteLock()
lib_get_algorithm(id::Lib_ID)::Optional{MarkovAlgorithm} = lock_read(LIB_ALGORITHM_LOCKER) do
    return get(LIB_ALGORITHMS, id, nothing)
end
lib_add_algorithm(algo::MarkovAlgorithm)::Lib_ID = lock_write(LIB_ALGORITHM_LOCKER) do
    id = LIB_NEXT_ALGORITHM_ID[]
    LIB_ALGORITHMS[id] = algo

    LIB_NEXT_ALGORITHM_ID[] += one(Lib_ID)
    return id
end
lib_remove_algorithm(id::Lib_ID)::Bool = lock_write(LIB_ALGORITHM_LOCKER) do
    prev_count = length(LIB_ALGORITHMS)
    delete!(LIB_ALGORITHMS, id)
    return length(LIB_ALGORITHMS) == (prev_count - 1)
end

# Algorithm states are stored in the same way, for the same reason.
const LIB_ALGO_STATES = Dict{Lib_ID, MarkovAlgoState}()
const LIB_NEXT_ALGO_STATE_ID = Ref{Lib_ID}(1)
const LIB_ALGO_STATE_LOCKER = ReadWriteLock()
lib_get_algo_state(id::Lib_ID)::Optional{MarkovAlgoState} = lock_read(LIB_ALGO_STATE_LOCKER) do
    return get(LIB_ALGO_STATES, id, nothing)
end
lib_add_algo_state(state::MarkovAlgoState)::Lib_ID = lock_write(LIB_ALGO_STATE_LOCKER) do
    id = LIB_NEXT_ALGO_STATE_ID[]
    LIB_ALGO_STATES[id] = state

    LIB_NEXT_ALGO_STATE_ID[] += one(Lib_ID)
    return id
end
function lib_remove_algo_state(algo_id::Lib_ID, state_id::Lib_ID)::Bool
    algo = lib_get_algorithm(algo_id)
    if isnothing(algo)
        return false
    end

    lock_write(LIB_ALGO_STATE_LOCKER) do
        if haskey(LIB_ALGO_STATES, state_id)
            close(LIB_ALGO_STATES[state_id], algo)
            delete!(LIB_ALGO_STATES, state_id)
            return true
        else
            return false
        end
    end
end

LIB_SUPPRESS_STDERR::Bool = false

"Writes an error string according to our protocol, and returns the byte size (or 0 if the error went to stderr)"
function lib_error_handler(e, out_c_str::Ptr{Cchar}, c_str_buffer_size::Integer, prefix_if_stderr::String)::Int
    if out_c_str == C_NULL
        if !LIB_SUPPRESS_STDERR
            print(stderr, "JMarkovJunior: ", prefix_if_stderr, " ")
            showerror(stderr, e)
            println(stderr)
        end
        return 0
    else
        err_str = sprint(showerror, e)

        # Edge-case: user gave us basically no memory to write the error message.
        if c_str_buffer_size < 4
            if c_str_buffer_size > 0
                unsafe_store!(out_c_str, Cchar('\0'), 1)
                return 1
            else
                return 0
            end
        # If the error is too long, cut it off with an ellipsis.
        elseif length(err_str) > c_str_buffer_size - 1
            for (c, i) in zip(err_str, 1:(c_str_buffer_size - 4))
                unsafe_store!(out_c_str, Cchar(c), i)
            end
            unsafe_store!(out_c_str, Cchar('.'),  c_str_buffer_size - 3)
            unsafe_store!(out_c_str, Cchar('.'),  c_str_buffer_size - 2)
            unsafe_store!(out_c_str, Cchar('.'),  c_str_buffer_size - 1)
            unsafe_store!(out_c_str, Cchar('\0'), c_str_buffer_size - 0)
            return c_str_buffer_size
        # Otherwise, write the full error.
        else
            for (i, c) in enumerate(err_str)
                unsafe_store!(out_c_str, Cchar(c), i)
            end
            unsafe_store!(out_c_str, Cchar('\0'), length(err_str) + 1)
            return length(err_str) + 1
        end
    end
end

Base.@ccallable jmj_get_suppess_stderr()::Cint = LIB_SUPPRESS_STDERR 
Base.@ccallable jmj_set_suppess_stderr(b::Cint)::Cvoid = (global LIB_SUPPRESS_STDERR = (b != 0))

Base.@ccallable function jmj_algo_parse(c_str::Cstring, out_err_str::Ptr{Cchar}, out_err_buffer_size::Cint)::Clonglong
    julia_str = unsafe_string(c_str)
    try
        return convert(Clonglong, lib_add_algorithm(markov_algo_parse(julia_str)))
    catch e
        return -convert(Clonglong,
            lib_error_handler(e, out_err_str, out_err_buffer_size,
                              "Failed to parse algorithm!")
        )
    end
end
Base.@ccallable function jmj_algo_close(id::Lib_ID)::Cint
    if lib_remove_algorithm(id)
        return 1
    else
        LIB_SUPPRESS_STDERR || println(stderr,
            "JMarkovJunior: No algorithm existed with the ID ", id, "! Nothing changes"
        )
        return 0
    end
end

Base.@ccallable function jmj_start(algo_id::Lib_ID,
                                   n_dims::Cint, size::Ptr{Cint},
                                   n_seed_bytes::Cint, seed_bytes::Ptr{UInt8})::Lib_ID
    # Validate and process size params.
    if n_dims < 1
        LIB_SUPPRESS_STDERR || println(stderr, "JMarkovJunior: Dimension must be at least 1; got ", n_dims)
        return zero(Lib_ID)
    end
    size_tuple = ntuple(i -> unsafe_load(size, i), n_dims)
    if any(size_tuple .< 1)
        LIB_SUPPRESS_STDERR || println(stderr, "JMarkovJunior: Size must be at least 1 along every axis; got ", size_tuple)
        return zero(Lib_ID)
    end

    # Validate and retrieve the algorithm.
    algo = lib_get_algorithm(algo_id)
    if isnothing(algo)
        LIB_SUPPRESS_STDERR || println(stderr, "JMarkovJunior: Unable to find MarkovJunior algorithm with ID ", algo_id, "!")
        return zero(Lib_ID)
    end

    # Retrive the seeds.
    # It's easier on the JIT if we don't have a ton of seed elements going into the PRNG,
    #    so prefer reinterpreting them as large uints.
    seeds = if n_seed_bytes < 0
        LIB_SUPPRESS_STDERR || println(stderr, "JMarkovJunior: Passed a negative number of RNG seeds (", n_seed_bytes, ")!")
        return zero(Lib_ID)
    elseif n_seed_bytes == 0
        rand(UInt64)
    elseif (n_seed_bytes % 8) == 0
        map(i -> unsafe_load(Ptr{UInt64}(seed_bytes), i),
            1:(n_seed_bytes ÷ 8))
    elseif (n_seed_bytes % 4) == 0
        map(i -> unsafe_load(Ptr{UInt32}(seed_bytes), i),
            1:(n_seed_bytes ÷ 4))
    elseif (n_seed_bytes ÷ 2) == 0
        map(i -> unsafe_load(Ptr{UInt16}(seed_bytes), i),
            1:(n_seed_bytes ÷ 2))
    else
        map(i -> unsafe_load(seed_bytes, i),
            1:n_seed_bytes)
    end

    algo_state = try
        markov_algo_start(algo, size_tuple, seeds)
    catch e
        if !LIB_SUPPRESS_STDERR
            print(stderr, "JMarkovJunior: Failed to start algorithm! ")
            showerror(stderr, e)
            println(stderr)
        end
        return zero(Lib_ID)
    end

    algo_state_id = lib_add_algo_state(algo_state)
    return algo_state_id
end
Base.Base.@ccallable function jmj_destroy(algo_id::Lib_ID, state_id::Lib_ID)::Cint
    if lib_remove_algo_state(algo_id, state_id)
        return 1
    else
        LIB_SUPPRESS_STDERR || println(stderr,
            "JMarkovJunior: No algo/state existed with the IDs",
            algo_id, "/", state_id, "! Nothing changes"
        )
        return 0
    end
end

Base.@ccallable function jmj_step(algo_id::Lib_ID, state_id::Lib_ID, count::Cint,
                                  out_was_error::Ptr{Cint})::Cint
    algo = lib_get_algorithm(algo_id)
    if isnothing(algo)
        LIB_SUPPRESS_STDERR || println(stderr, "JMarkovJunior: Unable to find MarkovJunior algorithm with ID ", algo_id, "!")
        (out_was_error == C_NULL) || unsafe_store!(out_was_error, one(Cint))
        return false
    end

    state = lib_get_algo_state(state_id)
    if isnothing(state)
        LIB_SUPPRESS_STDERR || println(stderr, "JMarkovJunior: Unable to find MarkovJunior state with ID ", state_id, "!")
        (out_was_error == C_NULL) || unsafe_store!(out_was_error, one(Cint))
        return false
    end

    # Hide the internal "-1 means forever" logic from the user.
    count = max(count, zero(Cint))

    # Stabilize the state type for the rest of the function.
    function finish(state::TState) where {TState<:MarkovAlgoState}
        markov_algo_step(algo, state, count)
        (out_was_error == C_NULL) || unsafe_store!(out_was_error, zero(Cint))
        return markov_algo_is_finished(algo, state)
    end
    return finish(state)
end
Base.@ccallable function jmj_finish(algo_id::Lib_ID, state_id::Lib_ID, out_was_error::Ptr{Cint})::Cvoid
    algo = lib_get_algorithm(algo_id)
    if isnothing(algo)
        LIB_SUPPRESS_STDERR || println(stderr, "JMarkovJunior: Unable to find MarkovJunior algorithm with ID ", algo_id, "!")
        (out_was_error == C_NULL) || unsafe_store!(out_was_error, one(Cint))
        return nothing
    end

    state = lib_get_algo_state(state_id)
    if isnothing(state)
        LIB_SUPPRESS_STDERR || println(stderr, "JMarkovJunior: Unable to find MarkovJunior state with ID ", state_id, "!")
        (out_was_error == C_NULL) || unsafe_store!(out_was_error, one(Cint))
        return nothing
    end

    (out_was_error == C_NULL) || unsafe_store!(out_was_error, zero(Cint))
    markov_algo_finish(algo, state)
    return nothing
end
Base.@ccallable function jmj_is_finished(algo_id::Lib_ID, state_id::Lib_ID, out_was_error::Ptr{Cint})::Cint
    algo = lib_get_algorithm(algo_id)
    if isnothing(algo)
        LIB_SUPPRESS_STDERR || println(stderr, "JMarkovJunior: Unable to find MarkovJunior algorithm with ID ", algo_id, "!")
        (out_was_error == C_NULL) || unsafe_store!(out_was_error, one(Cint))
        return false
    end

    state = lib_get_algo_state(state_id)
    if isnothing(state)
        LIB_SUPPRESS_STDERR || println(stderr, "JMarkovJunior: Unable to find MarkovJunior state with ID ", state_id, "!")
        (out_was_error == C_NULL) || unsafe_store!(out_was_error, one(Cint))
        return false
    end

    (out_was_error == C_NULL) || unsafe_store!(out_was_error, zero(Cint))
    return markov_algo_is_finished(algo, state)
end

Base.@ccallable function jmj_grid(state_id::Lib_ID,
                                  out_ndims::Ptr{Cint},
                                  out_size::Ptr{Cint}, size_capacity::Int)::Ptr{Cuchar}
    state = lib_get_algo_state(state_id)
    if isnothing(state)
        LIB_SUPPRESS_STDERR || println(stderr, "JMarkovJunior: Unable to find MarkovJunior state with ID ", state_id, "!")
        return C_NULL
    end

    # Stabilize the state type for the rest of the function.
    function finish(state::TState)::Ptr{UInt8} where {TState<:MarkovAlgoState}
        grid = markov_algo_grid(state)

        (out_ndims == C_NULL) || unsafe_store!(out_ndims, convert(Cint, ndims(grid)))
        (out_size == C_NULL) || for i in 1:min(ndims(grid), size_capacity)
            unsafe_store!(out_size, convert(Cint, size(grid, i)), i)
        end
        return pointer(grid)
    end
    return finish(state)
end
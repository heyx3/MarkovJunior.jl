# Start the IPC service and then try to use it.
# This doubles as a test of the DLL API.

#TODO: @suppress_err and then somehow log it when catching exceptions

(() -> begin # Wrap the test in a function to avoid global weirdness

PIPE_PATH = if Sys.iswindows()
    "\\\\.\\pipe\\jmj_ipc_unittest"
else
    "/tmp/jmj_ipc_unittest"
end
server = markovjunior_run_ipc(
    false, Val(false),
    allow_kill_message = true,
    pipe_path = PIPE_PATH
)
sleep(5) #TODO: Some kind of signal to know that the server has started
println("  (if you see this before the server says it's ready to accept clients, the test will fail due to race condition)")
channel = connect(PIPE_PATH)

try # Even if below tests fail, make sure to close the client and kill the server

# Do the handshake.
my_name = codeunits("Julia unit test -$(rand(1:1000))-")
write(channel, convert(UInt32, length(my_name)))
write(channel, my_name)

# Define helpers to call the IPC functions.
function ipc_parse_algo(algo_str::String, expected_algo_id::Integer, expect_success::Bool)::UInt32
    write(channel, UInt32(1))

    algo_bytes = codeunits(algo_str)
    write(channel, convert(UInt32, length(algo_bytes)))
    write(channel, algo_bytes)

    err_code = read(channel, UInt8)
    if err_code == 0
        err_str_len = read(channel, UInt32)
        err_str_bytes = read(channel, err_str_len)
        @bp_check(iszero(err_str_bytes[end]), "Error msg string didn't end in null-terminator")
        pop!(err_str_bytes)
        err_str = String(err_str_bytes)

        @bp_check(!expect_success,
                  "Failed to parse algorithm: (", err_code, ") ",
                    err_str, "\n", algo_str, "\n")
        println(stderr, "VVVVVVVVVVVVVVVVVVVVVVV\n",
                        "Error message from successful test, for reference: ", err_str,
                        "\n^^^^^^^^^^^^^^^^^^^^^")
        return zero(UInt32)
    elseif (err_code == 1)
        algo_id::UInt32 = read(channel, UInt32)
        if algo_id != expected_algo_id
            error("Expected parsed algo's ID to be ",
                  convert(Int, expected_algo_id), " but got ", convert(Int, algo_id))
        end
        @bp_check(expect_success, "Expected algo parse to fail but it did not! ", algo_str, "\n")
        return algo_id
    else
        error("Unexpected error code: ", err_code)
    end
end
function ipc_delete_algo(id::Integer, expect_success::Bool)::Nothing
    write(channel, UInt32(2))
    write(channel, convert(UInt32, id))
    err_code = read(channel, UInt8)
    if err_code == 1
        if expect_success
            # No problem!
        else
            error("Expected deleting algo ", id, " to fail but it worked somehow")
        end
    elseif err_code == 0
        if expect_success
            error("Failed to delete algo ", id)
        else
            # No problem!
        end
    else
        error("Unexpected error code: ", err_code)
    end
    return nothing
end
function ipc_start(algo_id::Integer, grid_size::Tuple{Vararg{Integer}}, seeds,
                   expect_size_allowed::Bool,
                   initial_state::Optional{Array},
                   expected_state_id::Integer, expect_success::Bool
                  )::UInt32
    grid_size_u32 = convert.(Ref(UInt32), grid_size)
    seed_bytes::Vector{UInt8} = reinterpret_bytes_slow(seeds)
    initial_state = if isnothing(initial_state)
        nothing
    else
        reshape(convert.(Ref(UInt8), initial_state),
                grid_size)
    end

    write(channel, UInt32(3))
    write(channel, convert(UInt32, algo_id))
    write(channel, convert(UInt32, length(grid_size)))
    write(channel, grid_size_u32...)

    err_code = read(channel, UInt8)
    if err_code == 1
        @bp_check(expect_size_allowed, "Grid size ", grid_size, " not allowed")

        if exists(initial_state)
            write(channel, one(UInt8))
            write(channel, initial_state)
        else
            write(channel, zero(UInt8))
        end

        write(channel, convert(UInt32, length(seed_bytes)))
        write(channel, seed_bytes)

        err_code = read(channel, UInt8)
        if err_code == 1
            @bp_check(expect_success, "ipc_start succeeded when it shouldn't have!")
            state_id = read(channel, UInt32)
            @bp_check(state_id == expected_state_id,
                    "Expected state ID to be ", expected_state_id, " but got ", state_id)
            return state_id
        elseif err_code == 0
            @bp_check(!expect_success, "ipc_start failed!")
            return zero(UInt32)
        else
            error("Unexpected error code: ", err_code)
        end
    else
        @bp_check(!expect_size_allowed, "Grid size ", grid_size, " should have been forbidden")
        return zero(UInt32)
    end
end
function ipc_destroy(state_id::Integer, expect_success::Bool)::Nothing
    write(channel, UInt32(4))
    write(channel, convert(UInt32, state_id))

    err_code = read(channel, UInt8)
    if err_code == 1
        @bp_check(expect_success, "ipc_destroy succeeded when it shouldn't have!")
    elseif err_code == 0
        @bp_check(!expect_success, "ipc_destroy failed!")
    else
        error("Unexpected error code: ", err_code)
    end
    return nothing
end
function ipc_advance(state_id::Integer,
                     # First int is 'lowest priority tick to count', second is '# ticks'
                     mode::Union{Val{:tags}, Val{:completed}, NTuple{2, Integer}},
                     expect_success::Bool, expect_finished::Bool,
                     expect_tag::Optional{Symbol}
                    )::Nothing
    write(channel, UInt32(5))
    write(channel, convert(UInt32, state_id))

    if mode isa NTuple{2, Integer}
        write(channel, UInt8(0))
        write(channel, convert(UInt32, mode[1]))
        write(channel, convert(UInt32, mode[2]))
    elseif tag_mode isa Val{:tags}
        write(channel, UInt8(1))
    elseif tag_mode isa Val{:completed}
        write(channel, UInt8(2))
    else
        error("Unhandled test case: ", mode)
    end

    err_code = read(channel, UInt8)
    if err_code == 0
        @bp_check(!expect_success, "ipc_advance() failed! Possibly invalid state ID ", state_id)
        return nothing
    elseif err_code != 1
        error("Unexpected error code: ", err_code)
    end
    @bp_check(expect_success, "ipc_advance() should have failed but it didn't! State ID ", state_id)

    algo_finished = read(channel, UInt8)
    if algo_finished == 0
        @bp_check(!expect_finished, "The algorithm should have finished here!")
        return nothing
    elseif algo_finished != 1
        error("Unexpected 'is algo finished' code: ", algo_finished)
    end
    @bp_check(expect_finished, "The algorithm shouldn't have finished here but it did!")

    had_tagged_event = read(channel, UInt8)
    if had_tagged_event == 0
        @bp_check(isnothing(expect_tag), "The algorithm should have hit a tagged event!")
        return nothing
    elseif had_tagged_event != 1
        error("Unexpected 'had tagged event' code: ", had_tagged_event)
    end

    tagged_event_bytes = Vector{UInt8}(undef, read(channel, UInt32))
    read!(channel, tagged_event_bytes)
    tagged_event = Symbol(tagged_event_bytes)
    @bp_check(exists(expect_tag), "The algorithm should not have hit a tagged event, but it hit '", tagged_event, "'")
    @bp_check(expect_tag == tagged_event, "Expected tagged event '", expect_tag, "' but got '", tagged_event, "'")

    return nothing
end
function ipc_grid(state_id::Integer, expect_success::Bool)::Optional{Array{UInt8}}
    write(channel, UInt32(6))
    write(channel, convert(UInt32, state_id))

    err_code = read(channel, UInt8)
    if err_code == 1
        @bp_check(expect_success, "ipc_grid succeeded when it shouldn't have!")

        n_dims = convert(Int, read(channel, UInt32))
        grid_size = map(i -> read(channel, UInt32), 1:n_dims)
        grid = Array{UInt8, n_dims}(undef, grid_size...)
        read!(channel, grid)

        return grid
    elseif err_code == 0
        @bp_check(!expect_success, "ipc_grid failed!")
        return nothing
    else
        error("Unexpected error code: ", err_code)
    end
end
function ipc_kill(expect_success::Bool)::Nothing
    write(channel, UInt32(7))

    err_code = read(channel, UInt8)
    if err_code == 1
        @bp_check(expect_success, "Killed the server even though we shouldn't have been able to!")
    elseif err_code == 0
        @bp_check(!expect_success, "Failed to kill the server!")
    else
        error("Unexpected error code: ", err_code)
    end

    return nothing
end

# Parse a new algorithm.
# Sprinkle in some bogus delete ops that should fail.
# Eventually delete the parsed algo for real.
ipc_delete_algo(1, false)
ipc_parse_algo("""@markovjunior 2 begin
    @rewrite 30 b=>w
end""", 1, true)
ipc_delete_algo(1, true)
ipc_delete_algo(1, false)
ipc_delete_algo(2, false)

# Send a bogus algorithm, then parse a real one for us to really test with.
ipc_parse_algo("""@markovjunior 2 begin
    @not_a_real_Op 42
end""", 0, false)
ipc_parse_algo("""@markovjunior 2 begin
    @rewrite 30 b=>w
end""", 2, true)
# Deleting with the old algo's ID should still fail.
ipc_delete_algo(1, false)

# At this point we should be able to kill the server, preventing new connections
#    while still using our own client connection.
ipc_kill(true)
sleep(1)
@bp_check(!isopen(server), "Server wasn't killed after 1 second")

# Try starting the algorithm.
# Also try some false starts.
ipc_start(1, (6, 6), (1, 4.5),   true, nothing,   0, false) # Failed due to algo ID
ipc_start(2, (4, ), (1, 4.5),  true, nothing,  0, false) # Failed due to algo being 2D and grid being 1D
ipc_start(2, ntuple(i -> Int(ceil(sqrt(MJ.IPC_DEFAULT_MAX_GRID_BYTE_SIZE)) + 1), 2),
          (1, 4.5),
          false, nothing,
          0, false) # Failed due to memory cap
ipc_start(2, (3, 12), (1, 4.5),   true, nothing,    1, true)
# Get the grid for the first time, and verify it.
let g = ipc_grid(1, true)
    @bp_check(size(g) == (3, 12), "Grid is ", size(g))
    @bp_check(all(iszero, g), "Grid: ", g)
end
ipc_grid(2, false)

# Run some iterations and check that there are now changed pixels.
ipc_advance(1, (3, 1), true, false, nothing)
let grid = ipc_grid(1, true)
    @bp_check(size(grid) == (3, 12), "Grid is ", size(grid))
    @bp_check(count(i->i==2, grid) == 1, "Grid: ", grid)
    @bp_check(count(iszero, grid) == 35, "Grid: ", grid)
end
ipc_advance(2, (3, 3), false, false, nothing) # Failed due to state ID
ipc_advance(1, (3, 3), true, false, nothing)
let grid = ipc_grid(1, true)
    @bp_check(size(grid) == (3, 12), "Grid is ", size(grid))
    @bp_check(count(i->i==2, grid) == 4, "Grid: ", grid)
    @bp_check(count(iszero, grid) == 32, "Grid: ", grid)
end
ipc_advance(1, Val(:completed), true, true, nothing)
let grid = ipc_grid(1, true)
    @bp_check(size(grid) == (3, 12), "Grid is ", size(grid))
    @bp_check(count(i->i==2, grid) == 30, "Grid: ", grid)
    @bp_check(count(iszero, grid) == 6, "Grid: ", grid)
end

# Destroy the state (and check related failure modes).
ipc_destroy(2, false)
ipc_destroy(1, true)
ipc_destroy(1, false)
ipc_grid(1, false)

# Verify the "finish" message and ability to write an initial state.
ipc_start(2, (4, 2), (1, 4.5),  true, [ 0, 1, 1, 0, 1, 0, 1, 0 ],    2, true)
let grid = ipc_grid(2, true)
    @bp_check(size(grid) == (4, 2), "Grid is ", size(grid))
    @bp_check(grid == UInt8[
        0 1
        1 0
        1 1
        0 0
    ], "Grid: ", grid)
end
ipc_advance(1, Val(:tags), false, false, nothing) # Failed due to state ID
ipc_advance(2, Val(:tags), true, true, nothing) # Finishes due to no tagged events
ipc_grid(1, false)
let grid = ipc_grid(2, true)
    @bp_check(size(grid) == (4, 2), "Grid is ", size(grid))
    @bp_check(grid == UInt8[
        2 1
        1 2
        1 1
        2 2
    ], "Grid: ", grid)
end
ipc_destroy(2, true)
ipc_destroy(2, false)

close(channel)
catch e
    # Wait a bit so the server thread can print its own stuff without interference.
    sleep(2)
    close(channel)
    sleep(2)
    rethrow()
end
end)()
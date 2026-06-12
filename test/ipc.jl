# Start the IPC service and then try to use it.
# This doubles as a test of the DLL API.

#TODO: @suppress_err and then somehow log it when catching exceptions

(() -> begin # Wrap the test in a function to avoid global weirdness

server = markovjunior_run_ipc(false, allow_kill_message=true)
sleep(5) #TODO: Some kind of signal to know that the server has started
println("  (if you see this before the server says it's ready to accept clients, the test will fail due to race condition)")
channel = connect(MarkovJunior.IPC_PIPE_PATH)

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
                   expected_state_id::Integer, expect_success::Bool
                  )::UInt32
    grid_size_u32 = convert.(Ref(UInt32), grid_size)
    seed_bytes::Vector{UInt8} = reinterpret_bytes_slow(seeds)

    write(channel, UInt32(3))
    write(channel, convert(UInt32, algo_id))
    write(channel, convert(UInt32, length(grid_size)))
    write(channel, grid_size_u32...)
    write(channel, convert(UInt32, length(seed_bytes)))
    write(channel, seed_bytes)

    err_code = read(channel, UInt8)
    if err_code == 1
        @bp_check(expect_success, error("ipc_start succeeded when it shouldn't have!"))
        state_id = read(channel, UInt32)
        @bp_check(state_id == expected_state_id,
                  "Expected state ID to be ", expected_state_id, " but got ", state_id)
        return state_id
    elseif err_code == 0
        @bp_check(!expect_success, error("ipc_start failed!"))
        return zero(UInt32)
    else
        error("Unexpected error code: ", err_code)
    end
end
function ipc_destroy(algo_id::Integer, state_id::Integer, expect_success::Bool)::Nothing
    write(channel, UInt32(4))
    write(channel, convert(UInt32, algo_id))
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
function ipc_step(algo_id::Integer, state_id::Integer, count::Integer,
                  expect_success::Bool, expect_finished::Bool
                 )::Bool
    write(channel, UInt32(5))
    write(channel, convert(UInt32, algo_id))
    write(channel, convert(UInt32, state_id))
    write(channel, convert(UInt32, count))

    err_code = read(channel, UInt8)
    if err_code == 1
        @bp_check(expect_success, "ipc_step succeeded when it shouldn't have!")
        is_algo_finished_code = read(channel, UInt8)
        if is_algo_finished_code == 1
            @bp_check(expect_finished, "ipc_step finished the algorithm when it shouldn't have!")
            return true
        elseif is_algo_finished_code == 0
            @bp_check(!expect_finished, "ipc_step didn't finish the algorithm like it should have!")
            return false
        else
            error("Expected 'is_algo_finished' to be 0 or 1, but it was ", is_algo_finished_code)
        end
    elseif err_code == 0
        @bp_check(!expect_success, "ipc_step failed!")
        return false
    else
        error("Unexpected error code: ", err_code)
    end
end
function ipc_finish(algo_id::Integer, state_id::Integer, expect_success::Bool)::Nothing
    write(channel, UInt32(6))
    write(channel, convert(UInt32, algo_id))
    write(channel, convert(UInt32, state_id))

    err_code = read(channel, UInt8)
    if err_code == 1
        @bp_check(expect_success, "ipc_finish succeeded when it shouldn't have!")
    elseif err_code == 0
        @bp_check(!expect_success, "ipc_finish failed!")
    else
        error("Unexpected error code: ", err_code)
    end
    return nothing
end
function ipc_is_finished(algo_id::Integer, state_id::Integer,
                         expect_success::Bool, expected_output::Bool)::Bool
    write(channel, UInt32(7))
    write(channel, convert(UInt32, algo_id))
    write(channel, convert(UInt32, state_id))

    err_code = read(channel, UInt8)
    if err_code == 1
        @bp_check(expect_success, "ipc_is_finished succeeded when it shouldn't have!")
        result = read(channel, UInt8)
        if result == 1
            @bp_check(expected_output, "ipc_is_finished returned true when it shouldn't have!")
            return true
        elseif result == 0
            @bp_check(!expected_output, "ipc_is_finished returned false when it shouldn't have!")
            return false
        else
            error("Unexpected return value from ipc_is_finished: ", result)
        end
    elseif err_code == 0
        @bp_check(!expect_success, "ipc_is_finished failed!")
        return false
    else
        error("Unexpected error code: ", err_code)
    end
end
function ipc_grid(state_id::Integer, expect_success::Bool)::Optional{Array{UInt8}}
    write(channel, UInt32(8))
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
    write(channel, UInt32(9))

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
ipc_start(1, (6, 6), (1, 4.5),   0, false) # Failed due to algo ID
ipc_start(2, ntuple(identity, 9), (1, 4.5),   0, false) # Failed due to 8D cap
ipc_start(2, (4, ), (1, 4.5),  0, false) # Failed due to algo being 2D and grid being 1D
ipc_start(2, ntuple(i -> 99999, 7), (1, 4.5),  0, false) # Failed due to memory cap)
ipc_start(2, (3, 12), (1, 4.5),    1, true)
# Get the grid for the first time, and verify it.
let g = ipc_grid(1, true)
    @bp_check(size(g) == (3, 12), "Grid is ", size(g))
    @bp_check(all(iszero, g), "Grid: ", g)
end
ipc_grid(2, false)

# Run some iterations and check that there are now changed pixels.
ipc_is_finished(2, 1, true, false)
ipc_step(2, 1, 1, true, false)
let grid = ipc_grid(1, true)
    @bp_check(size(grid) == (3, 12), "Grid is ", size(grid))
    @bp_check(count(i->i==2, grid) == 1, "Grid: ", grid)
    @bp_check(count(iszero, grid) == 35, "Grid: ", grid)
end
ipc_is_finished(2, 1, true, false)
ipc_step(1, 1, 1, false, false) # Failed due to algo ID
ipc_step(2, 2, 1, false, false) # Failed due to state ID
ipc_step(2, 1, 3, true, false)
let grid = ipc_grid(1, true)
    @bp_check(size(grid) == (3, 12), "Grid is ", size(grid))
    @bp_check(count(i->i==2, grid) == 4, "Grid: ", grid)
    @bp_check(count(iszero, grid) == 32, "Grid: ", grid)
end
ipc_step(2, 1, 100, true, true)
ipc_is_finished(2, 1, true, true)
ipc_is_finished(1, 1, false, false) # Failed due to algo ID
ipc_is_finished(2, 2, false, false) # Failed due to state ID
let grid = ipc_grid(1, true)
    @bp_check(size(grid) == (3, 12), "Grid is ", size(grid))
    @bp_check(count(i->i==2, grid) == 30, "Grid: ", grid)
    @bp_check(count(iszero, grid) == 6, "Grid: ", grid)
end

# Destroy the state (and check related failure modes).
ipc_destroy(2, 2, false)
ipc_destroy(1, 2, false)
ipc_destroy(2, 1, true)
ipc_destroy(2, 1, false)
ipc_grid(1, false)

# Verify the "finish" message.
ipc_start(2, (4, 2), (1, 4.5),    2, true)
ipc_finish(1, 2, false) # Failed due to algo ID
ipc_finish(2, 1, false) # Failed due to state ID
ipc_finish(2, 2, true)
ipc_is_finished(2, 2, true, true)
ipc_grid(1, false)
let grid = ipc_grid(2, true)
    @bp_check(size(grid) == (4, 2), "Grid is ", size(grid))
    @bp_check(count(i->i==2, grid) == 8, "Grid: ", grid)
end
ipc_destroy(2, 2, true)
ipc_destroy(2, 2, false)

finally
# End the server if it hasn't already.
if isopen(channel)
    write(channel, UInt32(9))
    @bp_check(read(channel, UInt8) in 0:1, "Final shutdown returned garbage?")
end
# End the client connection cleanly.
close(channel)
end
end)()
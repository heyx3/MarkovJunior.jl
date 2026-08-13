# Exposes RPC calls through a named pipe.
# Everything here is prefixed "ipc_".

"
The number at the end indicates the version;
  if we make a breaking change to the protocol then it will be incremented.
"
const IPC_PIPE_NAME = "jmj_ipc_2"
const IPC_PIPE_PATH = if Sys.iswindows()
    "\\\\.\\pipe\\$IPC_PIPE_NAME"
else
    "/tmp/$IPC_PIPE_NAME"
end

const IPC_Handle = UInt32

# Parsed algorithms are stored in a global lookup by unique ID.
const IPC_ALGO_LOOKUP = Dict{IPC_Handle, MarkovAlgorithm}()
const IPC_ALGORITHM_LOCKER = ReadWriteLock()
IPC_NEXT_ALGO_ID::IPC_Handle = 1
ipc_get_algorithm(id::IPC_Handle)::Optional{MarkovAlgorithm} = lock_read(IPC_ALGORITHM_LOCKER) do
    return get(IPC_ALGO_LOOKUP, id, nothing)
end
ipc_add_algorithm(algo::MarkovAlgorithm)::IPC_Handle = lock_write(IPC_ALGORITHM_LOCKER) do
    id = IPC_NEXT_ALGO_ID
    IPC_ALGO_LOOKUP[id] = algo

    IPC_NEXT_ALGO_ID += one(IPC_Handle)
    return id
end
ipc_remove_algorithm(id::IPC_Handle)::Optional{MarkovAlgorithm} = lock_write(IPC_ALGORITHM_LOCKER) do
    result = get(IPC_ALGO_LOOKUP, id, nothing)
    delete!(IPC_ALGO_LOOKUP, id)
    return result
end

# Running algorithms are managed in a similar way.
const IPC_STATE_LOOKUP = Dict{IPC_Handle, Tuple{MarkovAlgorithm, CellGrid, AlgoCommsChannel}}()
const IPC_STATE_LOCKER = ReadWriteLock()
IPC_NEXT_STATE_ID::IPC_Handle = 1
ipc_get_state(id::IPC_Handle)::Optional{Tuple{MarkovAlgorithm, CellGrid, AlgoCommsChannel}} = lock_read(IPC_STATE_LOCKER) do
    return get(IPC_STATE_LOOKUP, id, nothing)
end
ipc_add_state(algo::MarkovAlgorithm, @nospecialize(grid::CellGrid), channel::AlgoCommsChannel)::IPC_Handle = lock_write(IPC_STATE_LOCKER) do
    id = IPC_NEXT_STATE_ID
    IPC_STATE_LOOKUP[id] = (algo, grid, channel)

    IPC_NEXT_STATE_ID += one(IPC_Handle)
    return id
end
ipc_update_state_grid(id::IPC_Handle, @nospecialize(new_grid::CellGrid))::Bool = lock_write(IPC_STATE_LOCKER) do
    result = get(IPC_STATE_LOOKUP, id, nothing)
    if isnothing(result)
        return false
    end

    IPC_STATE_LOOKUP[id] = (result[0], new_grid, result[2])
    return true
end
ipc_remove_state(id::IPC_Handle)::Optional{Tuple{MarkovAlgorithm, CellGrid, AlgoCommsChannel}} = lock_write(IPC_ALGORITHM_LOCKER) do
    result = get(IPC_STATE_LOOKUP, id, nothing)
    delete!(IPC_STATE_LOOKUP, id)
    return result
end


"To prevent memory issues with garbage client data, we cap the max size of a grid"
const IPC_DEFAULT_MAX_GRID_BYTE_SIZE::Int = (1024 ^ 3) * 2
"To prevent memory issues with garbage client data, we cap the max size of a client's name"
const IPC_DEFAULT_MAX_CLIENT_NAME_BYTES::Int = 1024

const IPC_MAIN_START_CODE = UInt32(42)
const IPC_MAIN_STOP_CODE = UInt32(999)


"Internal macro to help with debug logging."
macro ipc_debug_log(args...)
    return :(
        $(esc(:DebugMode)) && println(stderr,
            $(esc(:client_name)), ": ",
            $(esc.(args)...)
        )
    )
end

const IPC_ERR_MSG_CUTOFF = UInt8[ codeunits("...")..., zero(UInt8) ]

function ipc_client_loop(client_name, channel, server, max_grid_byte_size, ::Val{DebugMode}) where {DebugMode}
    try
        while true
            msg_idx = read(channel, UInt32)
            # Message 1: Parse a new Algorithm
            if msg_idx == 1
                @ipc_debug_log client_name "|M: Parsing new Algorithm..."
                str_len = read(channel, UInt32)
                @ipc_debug_log client_name "|    algo is " str_len " bytes"
                str_bytes = read(channel, str_len)

                # Auto-append a null terminator if necessary.
                if !iszero(str_bytes[end])
                    @ipc_debug_log client_name "|    received with no null-terminator; appending it."
                    push!(str_bytes, zero(UInt8))
                end
                algo_src = unsafe_string(str_bytes)

                # Execute the parser.
                result = try
                    markov_algo_parse(algo_src)
                catch e
                    "Falied to parse algorithm! $(sprint(showerror, e))"
                end

                # Send/log the results.
                if result isa MarkovAlgorithm
                    @ipc_debug_log client_name "|        Parsing successful"
                    write(channel, one(UInt8))

                    handle = ipc_add_algorithm(result)
                    write(channel, handle)
                else
                    @ipc_debug_log client_name "|        " result
                    write(channel, zero(UInt8))

                    # Write the error message with a null-terminator.
                    err_str_bytes = codeunits(result)
                    write(channel, convert(UInt32, length(err_str_bytes) + 1))
                    write(channel, err_str_bytes)
                    write(channel, zero(UInt8))
                end
            # Message 2: Close a parsed Algorithm (not affecting current runs)
            elseif msg_idx == 2
                @ipc_debug_log client_name "|M: Closing a parsed Algorithm..."
                algo_id = read(channel, IPC_Handle)
                @ipc_debug_log client_name "|    Algo " algo_id

                result = ipc_remove_algorithm(algo_id)
                @ipc_debug_log client_name "|    result: " exists(result)
                write(channel, convert(UInt8, exists(result)))
            # Message 3: Start an algorithm run
            elseif msg_idx == 3
                @ipc_debug_log client_name "|M: Start running an Algorithm..."

                algo_id = read(channel, IPC_Handle)
                @ipc_debug_log client_name "|    Algo " algo_id

                n_dims = read(channel, UInt32)
                a_size = Vector{UInt32}(undef, n_dims)
                @ipc_debug_log client_name "|    Dims " n_dims
                read!(channel, a_size)
                @ipc_debug_log client_name "|    Size " Tuple(a_size)
                a_size = reinterpret(Cint, a_size)

                if prod(convert.(Ref(Int), a_size)) > max_grid_byte_size
                    println(stderr, client_name, "|    ERROR: requested too large a grid (", a_size, " == ", prod(a_size), ") and the request will fail")
                    write(channel, zero(UInt8))
                else
                    write(channel, one(UInt8))

                    # Read initial grid state (if applicable).
                    has_initial_grid = !iszero(read(channel, UInt8))
                    initial_grid::Optional{Array{UInt8, convert(Int, n_dims)}} = if has_initial_grid
                        @ipc_debug_log client_name "|    Provided an initial grid; reading now"
                        #TODO: Re-use an allocation for this
                        a = Array{UInt8, convert(Int, n_dims)}(undef, a_size...)
                        read!(channel, a)
                        a
                    else
                        nothing
                    end

                    # Read RNG seeds.
                    n_seed_bytes = read(channel, UInt32)
                    @ipc_debug_log client_name "|    Seed bytes " n_seed_bytes
                    seed_bytes::Vector{UInt8} = if n_seed_bytes == 0
                        [ rand(UInt8), rand(UInt8), rand(UInt8), rand(UInt8) ]
                    else
                        read(channel, n_seed_bytes)
                    end
                    @ipc_debug_log client_name "|    Seed: " Tuple(seed_bytes)

                    println(stderr, "#TODO: Read tick settings")
                    ticking = MarkovTickSettings()

                    # Execute the start operation.
                    algo = ipc_get_algorithm(algo_id)
                    if !exists(algo)
                        println(stderr, client_name, "|    ERROR: Algorithm handle ", algo_id, " not found")
                        write(channel, zero(UInt8))
                    else
                        try
                            run_channel::AlgoCommsChannel = markov_algo_run(
                                algo,
                                has_initial_grid ? initial_grid : Tuple(a_size),
                                ticking,
                                seeds=seed_bytes
                            )

                            # Move the algorithm forward until we've got a reference to the initial grid.
                            grid::Optional{CellGrid} = nothing
                            while isnothing(grid)
                                next_tick = take!(run_channel)
                                if next_tick == TAG_NEW_GRID
                                    grid = take!(run_channel)
                                end
                            end

                            handle = ipc_add_state(algo, grid, run_channel)

                            @ipc_debug_log client_name "|    Run handle is " handle
                            write(channel, one(UInt8))
                            write(channel, handle)
                        catch e
                            println(stderr, client_name, "|    ERROR: ", sprint(showerror, e))
                            write(channel, zero(UInt8))
                        end
                    end
                end
            # Message 4: Destroy a running state
            elseif msg_idx == 4
                @ipc_debug_log client_name "|M: Destroy a running state..."
                state_id = read(channel, IPC_Handle)
                @ipc_debug_log client_name "|     State=" state_id

                result = ipc_remove_state(state_id)
                if exists(result)
                    @ipc_debug_log client_name "|     Successful"
                    write(channel, one(UInt8))
                    (algo, grid, comms_channel) = result
                    close(comms_channel) # Will trigger an exception upon the algo's next tick,
                                         #    which is then handled appropriately
                else
                    println(stderr, client_name, "|    ERROR: state with that ID doesn't exist")
                    write(channel, zero(UInt8))
                end
            # Message 5: Advance the algorithm forward
            elseif msg_idx == 5
                @ipc_debug_log client_name "|M: Step forward..."
                state_id = read(channel, UInt32)
                @ipc_debug_log client_name "|     state=" state_id

                # Get the type of tick logic to run.
                tick_type_code = read(channel, UInt8)
                (tag_cutoff, tick_cutoff, tick_count)::Tuple{Symbol, Int, Int} = if tick_type_code == 0 # Normal tick
                    cutoff_priority = read(channel, UInt32)
                    tick_count = read(channel, UInt32)
                    @ipc_debug_log client_name "|    " tick_count " ticks of " cutoff_priority
                    (:all, cutoff_priority, tick_count)
                elseif tick_type_code == 1 # Tagged event
                    @ipc_debug_log client_name "|    tagged events"
                    (:all, typemax(Int), typemax(Int))
                elseif tick_type_code == 2 # Run to completion
                    @ipc_debug_log client_name "|    completion"
                    (:completion, typemax(Int), typemax(Int))
                else
                    println(stderr, client_name, "|    Invalid tick mode ", tick_type_code)
                    (:error, typemax(Int), typemax(Int))
                end

                # Execute the tick(s).
                state_query_result = ipc_get_state(state_id)
                if isnothing(state_query_result)
                    println(stderr, client_name, "|    Unknown state ID: ", state_id)
                    write(channel, zero(UInt8))
                elseif tag_cutof == :error
                    write(channel, zero(UInt8))
                else
                    (algo, grid, state_channel) = state_query_result
                    algo_finished::Bool = false
                    user_request_finished::Bool = false
                    reported_tag::Optional{Symbol} = nothing
                    while true
                        # Tick and accept the result.
                        put!(state_channel, zero(Int))
                        next_tick_data = take!(state_channel)

                        if next_tick_data isa Int
                            if next_tick_data >= tick_cutoff
                                tick_count -= 1
                                if tick_count <= 0
                                    @ipc_debug_log client_name "|    Hit tick limit!"
                                    user_request_finished = true
                                    break
                                end
                            end
                        elseif next_tick_data isa Symbol
                            # If the algorithm was just emitting a 'start' tag,
                            #    ignore it.
                            if next_tick_data != TAG_ALGO_STARTING
                                # If the algorithm completed, we must finish the tick.
                                if next_tick_data == TAG_ALGO_COMPLETED
                                    algo_finished = true
                                    user_request_finished = tag_cutoff in (:all, :completion)
                                    reported_tag = next_tick_data
                                    break
                                end
                                # If the algorithm reallocated the grid, store the new grid.
                                if next_tick_data == TAG_NEW_GRID
                                    grid = take!(state_channel)
                                    ipc_update_state_grid(state_id, grid)
                                end
                                # If the user wanted any tagged event, then we're done.
                                if tag_cutoff == :all
                                    user_request_finished = true
                                    reported_tag = next_tick_data
                                    break
                                end
                            end
                        else
                            error("Unhandled ", typeof(next_tick_data), ": ", next_tick_data)
                        end
                    end

                    # Report to the user what happened.
                    write(channel, one(UInt8))
                    write(channel, convert(UInt8, algo_finished))
                    if !algo_finished
                        write(channel, convert(UInt8, exists(reported_tag)))
                        if exists(reported_tag)
                            # Julia lets us write a Symbol's string representation directly to stream,
                            #    but without the null-terminator.
                            write(channel, convert(UInt32, sizeof(reported_tag) + 1))
                            write(channel, reported_tag)
                            write(channel, zero(UInt8))
                        end
                    end
                end

                if was_error == 1
                    @ipc_debug_log "    failed!"
                    write(channel, zero(UInt8))
                else
                    @ipc_debug_log "    result: " result
                    write(channel, one(UInt8))
                    write(channel, convert(UInt8, result))
                end
            # Message 6: Download the grid
            elseif msg_idx == 6
                @ipc_debug_log client_name "|M: Download the grid..."
                state_id = read(channel, UInt32)
                @ipc_debug_log "    State=" state_id

                result = ipc_get_state(state_id)
                if exists(result)
                    @ipc_debug_log client_name "|     Successful"
                    write(channel, one(UInt8))
                    (algo, grid, comms_channel) = result

                    write(channel, convert(UInt32, ndims(grid)))
                    for s in size(grid)
                        write(channel, convert(UInt32, s))
                    end
                    write(channel, grid)
                else
                    println(stderr, client_name, "|    ERROR: state with that ID doesn't exist")
                    write(channel, zero(UInt8))
                end
            # Message 7: Kill sever thread
            elseif msg_idx == 7
                @ipc_debug_log client_name "|M: Kill server thread..."
                if isnothing(server)
                    println(stderr, "Client \"", client_name, "\" asked to kill the server thread but that's not allowed.")
                    write(channel, zero(UInt8))
                elseif isopen(server)
                    println(stderr, "Client \"", client_name, "\" killed the server thread! No more clients can connect.")
                    close(server)
                    write(channel, one(UInt8))
                else
                    println(stderr, "Client \"", client_name, "\" tried to kill the server thread but it's already dead")
                    write(channel, one(UInt8))
                end
            else
                println(stderr, "Client \"", client_name, "\" sent invalid message index ", msg_idx,
                                "! Everything it does from now on is almost certainly garbage, ",
                                " so we're closing the connection.")
                break
            end
        end
    catch e
        if e isa EOFError
            println(stderr, "Client \"", client_name, "\" closed the connection.")
        else
            println(stderr, "Client \"", client_name,
                            "\" caused an exception! The communication protocol is broken so their connection will be cut off.")
            showerror(stderr, e, catch_backtrace())
            print(stderr, "\n\n")
        end
    finally
        if isopen(channel)
            close(channel)
        end
    end

    return nothing
end

"
Starts the IPC server.
If not blocking, the server object is returned
  which you can `close` at any time to prevent new clients.

The protocol is described in the project readme, puposefully not duplicated in our doc-strings.

It has the following optional keyword arguments:

* `allow_kill_message=true` for an extra IPC handler that allows you to kill the server thread (existing clients aren't affected).
* `max_grid_byte_size=x` to change the cap for a grid's memory allocation.
* `max_client_name_byte_size=x` to change the cap for a client's cosmetic name.
* `server_is_ready_callback = () -> ...` to react to the moment the server is ready to accept clients on the named pipe.
This happens relatively quickly and is always invoked on the calling thread, before this function returns.
"
function markovjunior_run_ipc(blocking::Bool, ::Val{DebugMode}
                              ;
                              max_grid_byte_size::Int=IPC_DEFAULT_MAX_GRID_BYTE_SIZE,
                              max_client_name_byte_size::Int=IPC_DEFAULT_MAX_CLIENT_NAME_BYTES,
                              allow_kill_message::Bool = false,
                              pipe_path::String = IPC_PIPE_PATH,
                              server_is_ready_callback = () -> println(stderr, "\tNow open to clients!")
                             ) where {DebugMode}
    println(stderr, "Offering JMarkovJunior at ", pipe_path)
    server = listen(pipe_path)
    server_is_ready_callback()

    server_loop = () -> while true
        # Note that clients can signal the entire server to end by simply closing it.
        channel = try
            accept(server)
        catch e
            if isopen(server)
                rethrow()
            else
                println(stderr, "Server noticed the kill signal; closing...")
                break
            end
        end

        # Initiate the handshake.
        client_name = try
            client_name_len = read(channel, UInt32)
            if client_name_len > max_client_name_byte_size
                println(stderr, "Client tried to send a name with ", client_name_len, " bytes in it!")
                error()
            else
                cn_array = read(channel, client_name_len)

                # Trim the null-terminator.
                if !isempty(cn_array) && cn_array[end] == 0
                    pop!(cn_array)
                end
                cn = String(cn_array)

                println(stderr, "New client: \"", cn, "\"")
                cn
            end
        catch e
            println(stderr, "New client failed their handshake and will be rejected.")
            close(channel)
            nothing
        end

        # Start the client's message loop in a new task.
        exists(client_name) && Threads.@spawn(ipc_client_loop(
            $client_name, $channel,
            $allow_kill_message ? server : nothing,
            max_grid_byte_size,
            Val(DebugMode)
        ))
    end
    if blocking
        server_loop()
        return nothing
    else
        Threads.@spawn try
            server_loop()
        catch e
            @error "IPC Server error!" exception=(e, catch_backtrace())
        end
        return server
    end
end

function markovjunior_run_ipc_main()::Cint
    arg_settings = ArgParseSettings()
    @add_arg_table! arg_settings begin
        "--max-grid", "-g"
            help = "The maximum number of bytes allowed for a single grid; default $(IPC_DEFAULT_MAX_GRID_BYTE_SIZE)"
            arg_type = Int
            default = IPC_DEFAULT_MAX_GRID_BYTE_SIZE
        "--max-client-name", "-c"
            help = "The maximum number of bytes allowed for a client's cosmetic name; default $(IPC_DEFAULT_MAX_CLIENT_NAME_BYTES)"
            arg_type = Int
            default = IPC_DEFAULT_MAX_CLIENT_NAME_BYTES
        "--allow-kill-msg", "-k"
            help = "If true, any client may send a special message to kill the server thread (preventing new clients)"
            action = :store_true
        "--timeout", "-t"
            help = "Server gives up and dies after this many seconds of zero client activity"
            arg_type = Float32
            default = +Inf
        "--debug", "-d"
            help = "Verbose logging of server logic to pinpoint crashes"
            action = :store_true
        "pipe_path"
            help = "An alternative name for the named-pipe that other processes talk to (valid syntax depends on your OS)"
            default = IPC_PIPE_PATH
    end
    arg_table = parse_args(ARGS, arg_settings)

    try
        markovjunior_run_ipc(
            true, Val(arg_table["debug"]::Bool),
            allow_kill_message = arg_table["allow-kill-msg"],
            max_grid_byte_size = arg_table["max-grid"],
            max_client_name_byte_size = arg_table["max-client-name"],
            pipe_path = arg_table["pipe_path"],
            server_is_ready_callback = () -> begin
                print(stderr, "\tNow open to clients! ")
                write(stdout, IPC_MAIN_START_CODE)
                flush(stdout)
                println(stderr)
            end
        )
        println(stderr)
        write(stdout, IPC_MAIN_STOP_CODE)
        flush(stdout)
        println(stderr, "Just wrote the stop code")
        println(stderr)
        return 0
    catch e
        showerror(stderr, e, catch_backtrace())
        write(stdout, IPC_MAIN_STOP_CODE)
        println(stderr, "Just wrote the stop code")
        return 1
    end
end
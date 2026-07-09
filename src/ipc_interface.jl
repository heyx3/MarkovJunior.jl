# Wraps the C-like API in a named pipe IPC setup.
# Everything is prefixed "ipc_".

"
The number at the end indicates the version;
  if we make a breaking change to the protocol then it will be incremented.
"
const IPC_PIPE_NAME = "jmj_ipc_1"
const IPC_PIPE_PATH = if Sys.iswindows()
    "\\\\.\\pipe\\$IPC_PIPE_NAME"
else
    "/tmp/$IPC_PIPE_NAME"
end

"For now we need to put a cap on the max dimensionality of the grid under the IPC protocol."
const IPC_MAX_DIMS = 8

"To prevent memory issues with garbage client data, we cap the max size of a grid"
const IPC_DEFAULT_MAX_GRID_BYTE_SIZE::Int = (1024 ^ 3) * 2
"To prevent memory issues with garbage client data, we cap the max size of a client's name"
const IPC_DEFAULT_MAX_CLIENT_NAME_BYTES::Int = 1024

const IPC_MAIN_START_CODE = UInt32(42)
const IPC_MAIN_STOP_CODE = UInt32(999)

#TODO: Client connections should remember the size of each grid so we can remove the 8D cap, and store the full set of their algorithms+states so it can auto-clean them up when the connection closes.


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
            if msg_idx == 1 # jmj_algo_parse
                @ipc_debug_log "jmj_algo_parse()..."
                str_len = read(channel, UInt32)
                @ipc_debug_log "    algo is " str_len " bytes"
                str_bytes = read(channel, str_len)
                # Auto-append a null terminator if necessary.
                if !iszero(str_bytes[end])
                    @ipc_debug_log "    received with no null-terminator; appending it."
                    push!(str_bytes, zero(UInt8))
                end

                err_str_buffer = Vector{Cchar}(undef, 1024)

                result = GC.@preserve str_bytes err_str_buffer begin
                    jmj_algo_parse(convert(Cstring, pointer(reinterpret(Cchar, str_bytes))),
                                   pointer(err_str_buffer), convert(Cint, length(err_str_buffer)))
                end
                @ipc_debug_log "    result: " result

                if result < 1
                    write(channel, zero(UInt8))
                    # Note that the error string includes a null terminator
                    #    because we called through the C API.
                    err_len = -result
                    @ipc_debug_log "        That means we have an error string of " err_len " bytes:"
                    @ipc_debug_log "            " String(reinterpret(UInt8, @view(err_str_buffer[1:err_len])))
                    write(channel, convert(UInt32, err_len))
                    write(channel, @view(err_str_buffer[1:err_len]))
                else
                    write(channel, one(UInt8))
                    write(channel, convert(UInt32, result))
                end
            elseif msg_idx == 2 # jmj_algo_close
                @ipc_debug_log "jmj_algo_close()..."
                algo_id = read(channel, UInt32)
                @ipc_debug_log "    Algo " algo_id
                result = jmj_algo_close(convert(Lib_ID, algo_id))
                @ipc_debug_log "    result: " result
                write(channel, convert(UInt8, result))
            elseif msg_idx == 3 # jmj_start
                @ipc_debug_log "jmj_start()..."

                algo_id = read(channel, UInt32)
                @ipc_debug_log "    Algo " algo_id

                n_dims = read(channel, UInt32)
                a_size = Vector{UInt32}(undef, n_dims)
                @ipc_debug_log "    Dims " n_dims
                read!(channel, a_size)
                @ipc_debug_log "    Size " Tuple(a_size)
                a_size = reinterpret(Cint, a_size)

                if n_dims > IPC_MAX_DIMS
                    println(stderr, "Client \"", client_name, "\" requested a ", n_dims, "D grid which is currently not allowed")
                    write(channel, zero(UInt8))
                elseif prod(convert.(Ref(Int), a_size)) > max_grid_byte_size
                    println(stderr, "Client \"", client_name, "\" requested too large a grid (", prod(a_size), ") and the request will fail")
                    write(channel, zero(UInt8))
                else
                    write(channel, one(UInt8))

                    initial_grid::Optional{Array{UInt8, convert(Int, n_dims)}} = if iszero(read(channel, UInt8))
                        nothing
                    else
                        #TODO: Re-use an allocation for this
                        a = Array{UInt8, convert(Int, n_dims)}(undef, a_size...)
                        read!(channel, a)
                        a
                    end

                    n_seed_bytes = read(channel, UInt32)
                    @ipc_debug_log "    Seed bytes " n_seed_bytes
                    seed_bytes = read(channel, n_seed_bytes)
                    @ipc_debug_log "    Seed: " Tuple(seed_bytes)

                    result = GC.@preserve a_size seed_bytes initial_grid begin
                        jmj_start(convert(Lib_ID, algo_id),
                                  convert(Cint, n_dims), pointer(a_size),
                                  convert(Cint, n_seed_bytes), pointer(seed_bytes),
                                  isnothing(initial_grid) ? Ptr{Cuchar}(C_NULL) : pointer(initial_grid))
                    end
                    if iszero(result)
                        write(channel, zero(UInt8))
                        @ipc_debug_log "    wrote failure code"
                    else
                        write(channel, one(UInt8))
                        @ipc_debug_log "    wrote success code; result is " result
                        write(channel, convert(UInt32, result))
                    end
                end
            elseif msg_idx == 4 # jmj_destroy
                @ipc_debug_log "jmj_destroy()..."
                algo_id = read(channel, UInt32)
                state_id = read(channel, UInt32)
                @ipc_debug_log "    Algo=" algo_id "   State=" state_id
                result = jmj_destroy(convert(Lib_ID, algo_id), convert(Lib_ID, state_id))
                @ipc_debug_log "    Result=" result
                write(channel, convert(UInt8, result))
            elseif msg_idx == 5 # jmj_step
                @ipc_debug_log "jmj_step()..."
                algo_id = read(channel, UInt32)
                @ipc_debug_log "    algo " algo_id
                state_id = read(channel, UInt32)
                @ipc_debug_log "    state " state_id
                count = read(channel, UInt32)
                @ipc_debug_log "    step count " count

                (result, was_error) = let was_error = Ref{Cint}()
                    r = GC.@preserve was_error begin
                        jmj_step(convert(Lib_ID, algo_id), convert(Lib_ID, state_id),
                                convert(Cint, count),
                                Base.unsafe_convert(Ptr{Cint}, was_error))
                    end
                    (r, was_error[])
                end

                if was_error == 1
                    @ipc_debug_log "    failed!"
                    write(channel, zero(UInt8))
                else
                    @ipc_debug_log "    result: " result
                    write(channel, one(UInt8))
                    write(channel, convert(UInt8, result))
                end
            elseif msg_idx == 6 # jmj_finish
                @ipc_debug_log "jmj_finish()..."
                algo_id = read(channel, UInt32)
                @ipc_debug_log "    Algo " algo_id
                state_id = read(channel, UInt32)
                @ipc_debug_log "    State " state_id

                results = @timed let was_error = Ref{Cint}()
                    GC.@preserve was_error begin
                        jmj_finish(convert(Lib_ID, algo_id), convert(Lib_ID, state_id),
                                Base.unsafe_convert(Ptr{Cint}, was_error))
                    end
                    was_error[]
                end
                was_error = results.value
                @ipc_debug_log "    " sprint(io -> Base.time_print(io,
                    # This is the internal Base code that prints @time.
                    results.time*1e9,
                    results.gcstats.allocd, results.gcstats.total_time, Base.gc_alloc_count(results.gcstats),
                    results.lock_conflicts,
                    results.compile_time*1e9, results.recompile_time*1e9,
                    false
                ))

                if was_error == 1
                    @ipc_debug_log "    failed!"
                    write(channel, zero(UInt8))
                else
                    @ipc_debug_log "    succeeded"
                    write(channel, one(UInt8))
                end
            elseif msg_idx == 7 # jmj_is_finished
                @ipc_debug_log "jmj_is_finished()..."
                algo_id = read(channel, UInt32)
                @ipc_debug_log "    Algo " algo_id
                state_id = read(channel, UInt32)
                @ipc_debug_log "    State " state_id

                (result, was_error) = let was_error = Ref{Cint}()
                    r = GC.@preserve was_error begin
                        jmj_is_finished(convert(Lib_ID, algo_id), convert(Lib_ID, state_id),
                                        Base.unsafe_convert(Ptr{Cint}, was_error))
                    end
                    (r, was_error[])
                end

                if was_error == 1
                    @ipc_debug_log "    failed!"
                    write(channel, zero(UInt8))
                else
                    @ipc_debug_log "    succeeded! Result: " result
                    write(channel, one(UInt8))
                    write(channel, convert(UInt8, result))
                end
            elseif msg_idx == 8 # jmj_grid
                @ipc_debug_log "jmj_grid()..."
                state_id = read(channel, UInt32)
                @ipc_debug_log "    State " state_id

                let dims_buf = Ref{Cint}(),
                    size_buf = Ref{NTuple{IPC_MAX_DIMS, Cint}}()
                  grid_ptr = GC.@preserve dims_buf size_buf begin
                        jmj_grid(convert(Lib_ID, state_id),
                                Base.unsafe_convert(Ptr{Cint}, dims_buf),
                                Base.unsafe_convert(Ptr{Cint}, size_buf),
                                IPC_MAX_DIMS)
                  end
                  if grid_ptr == C_NULL
                      @ipc_debug_log "    failed!"
                      write(channel, zero(UInt8))
                  else
                      @ipc_debug_log "    succeeded!"
                      write(channel, one(UInt8))

                      @ipc_debug_log "    Dims " dims_buf[]
                      write(channel, convert(UInt32, dims_buf[]))

                      @ipc_debug_log "    Size " size_buf[][1:dims_buf[]]
                      for i in 1:dims_buf[]
                          write(channel, convert(UInt32, size_buf[][i]))
                      end

                      n_grid_bytes = prod(size_buf[][i] for i in 1:dims_buf[])
                      @ipc_debug_log "   Grid bytes " n_grid_bytes
                      unsafe_write(channel, grid_ptr, n_grid_bytes)
                  end
                end
            elseif msg_idx == 9 # kill server thread
                if isnothing(server)
                    println(stderr, "Client \"", client_name, "\" asked to kill the server thread but that's not allowed.")
                    write(channel, zero(UInt8))
                elseif isopen(server)
                    println(stderr, "Client \"", client_name, "\" asked to kill the server thread! No more clients can connect.")
                    close(server)
                    write(channel, one(UInt8))
                else
                    println(stderr, "Client \"", client_name, "\" tried to kill the server but it's already dead")
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
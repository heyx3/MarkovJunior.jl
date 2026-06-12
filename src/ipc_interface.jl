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

#TODO: Add a timeout for finish operations to prevent the whole service from collapsing under load
#TODO: Make each specific client-connection more concrete. It should remember the size of each grid so we can remove the 8D cap, and store the full set of their algorithms+states so it can auto-clean them up when the connection closes.

"
An array of functions that clients can call in our IPC protocol.
The protocol is described in the project readme, puposefully not duplicated in our doc-strings.
"
const IPC_HANDLERS = tuple(
    # The "suppress stderr" flags aren't supported as there's no point -- this is its own process

    # 1: jmj_algo_parse
    (channel, client_name::String, max_grid_byte_size::Int) -> begin
        str_len = read(channel, UInt32)
        str_bytes = read(channel, str_len)
        # Auto-append a null terminator if necessary.
        if !iszero(str_bytes[end])
            push!(str_bytes, zero(UInt8))
        end

        err_str_buffer = Vector{Cchar}(undef, 1024)

        result = GC.@preserve str_bytes err_str_buffer begin
            jmj_algo_parse(convert(Cstring, pointer(str_bytes)),
                           pointer(err_str_buffer), convert(Cint, length(err_str_buffer)))
        end

        if result < 1
            write(channel, zero(UInt8))
            # Note that the error string includes a null terminator
            #    because we called through the C API.
            err_len = -result
            write(channel, convert(UInt32, err_len))
            write(channel, @view(err_str_buffer[1:err_len]))
        else
            write(channel, one(UInt8))
            write(channel, convert(UInt32, result))
        end

        return nothing
    end,
    # 2: jmj_algo_close
    (channel, client_name::String, max_grid_byte_size::Int) -> begin
        algo_id = read(channel, UInt32)
        result = jmj_algo_close(convert(Lib_ID, algo_id))
        write(channel, convert(UInt8, result))

        return nothing
    end,

    # 3: jmj_start
    (channel, client_name::String, max_grid_byte_size::Int) -> begin
        algo_id = read(channel, UInt32)

        n_dims = read(channel, UInt32)
        size = Vector{UInt32}(undef, n_dims)
        read!(channel, size)
        size = reinterpret(Cint, size)

        n_seed_bytes = read(channel, UInt32)
        seed_bytes = read(channel, n_seed_bytes)

        if n_dims > IPC_MAX_DIMS
            println(stderr, "Client \"", client_name, "\" requested a ", n_dims, "D grid which is currently not allowed")
            write(channel, zero(UInt8))
            return nothing
        elseif prod(convert.(Ref(Int), size)) > max_grid_byte_size
            println(stderr, "Client \"", client_name, "\" requested too large a grid (", prod(size), ") and the request will fail")
            write(channel, zero(UInt8))
            return nothing
        end

        result = GC.@preserve size seed_bytes begin
            jmj_start(convert(Lib_ID, algo_id),
                      convert(Cint, n_dims), pointer(size),
                      convert(Cint, n_seed_bytes), pointer(seed_bytes))
        end

        if iszero(result)
            write(channel, zero(UInt8))
        else
            write(channel, one(UInt8))
            write(channel, convert(UInt32, result))
        end

        return nothing
    end,
    # 4: jmj_destroy
    (channel, client_name::String, max_grid_byte_size::Int) -> begin
        algo_id = read(channel, UInt32)
        state_id = read(channel, UInt32)
        result = jmj_destroy(convert(Lib_ID, algo_id), convert(Lib_ID, state_id))
        write(channel, convert(UInt8, result))

        return nothing
    end,

    # 5: jmj_step
    (channel, client_name::String, max_grid_byte_size::Int) -> begin
        algo_id = read(channel, UInt32)
        state_id = read(channel, UInt32)
        count = read(channel, UInt32)
        (result, was_error) = let was_error = Ref{Cint}()
            r = GC.@preserve was_error begin
                jmj_step(convert(Lib_ID, algo_id), convert(Lib_ID, state_id),
                         convert(Cint, count),
                         Base.unsafe_convert(Ptr{Cint}, was_error))
            end
            (r, was_error[])
        end

        if was_error == 1
            write(channel, zero(UInt8))
        else
            write(channel, one(UInt8))
            write(channel, convert(UInt8, result))
        end

        return nothing
    end,
    # 6: jmj_finish
    (channel, client_name::String, max_grid_byte_size::Int) -> begin
        algo_id = read(channel, UInt32)
        state_id = read(channel, UInt32)
        was_error = let was_error = Ref{Cint}()
            GC.@preserve was_error begin
                jmj_finish(convert(Lib_ID, algo_id), convert(Lib_ID, state_id),
                           Base.unsafe_convert(Ptr{Cint}, was_error))
            end
            was_error[]
        end

        if was_error == 1
            write(channel, zero(UInt8))
        else
            write(channel, one(UInt8))
        end

        return nothing
    end,
    # 7: jmj_is_finished
    (channel, client_name::String, max_grid_byte_size::Int) -> begin
        algo_id = read(channel, UInt32)
        state_id = read(channel, UInt32)
        (result, was_error) = let was_error = Ref{Cint}()
            r = GC.@preserve was_error begin
                jmj_is_finished(convert(Lib_ID, algo_id), convert(Lib_ID, state_id),
                                Base.unsafe_convert(Ptr{Cint}, was_error))
            end
            (r, was_error[])
        end

        if was_error == 1
            write(channel, zero(UInt8))
        else
            write(channel, one(UInt8))
            write(channel, convert(UInt8, result))
        end

        return nothing
    end,

    # 8: jmj_grid
    (channel, client_name::String, max_grid_byte_size::Int) -> begin
        state_id = read(channel, UInt32)

        let dims_buf = Ref{Cint}(),
            size_buf = Ref{NTuple{IPC_MAX_DIMS, Cint}}()
          grid_ptr = GC.@preserve dims_buf size_buf begin
              jmj_grid(convert(Lib_ID, state_id),
                       Base.unsafe_convert(Ptr{Cint}, dims_buf),
                       Base.unsafe_convert(Ptr{Cint}, size_buf),
                       IPC_MAX_DIMS)
          end
          if grid_ptr == C_NULL
              write(channel, zero(UInt8))
          else
              write(channel, one(UInt8))
              write(channel, convert(UInt32, dims_buf[]))
              for i in 1:dims_buf[]
                  write(channel, convert(UInt32, size_buf[][i]))
              end
              unsafe_write(channel, grid_ptr, prod(size_buf[][i] for i in 1:dims_buf[]))
          end
        end
    end

    # 9: special "kill the server thread" message that has no official handler.
    # This call helps one program take complete ownership of this service's lifetime.
    # You can optionally disable this message when starting the server.
)

function ipc_client_loop(client_name, channel, server, max_grid_byte_size)
    try
        while true
            msg_idx = read(channel, UInt32)
            if msg_idx == length(IPC_HANDLERS) + 1
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
            elseif (msg_idx < 1) || (msg_idx > length(IPC_HANDLERS))
                println(stderr, "Client \"", client_name, "\" sent invalid message index ", msg_idx,
                                "! Everything it does from now on is almost certainly garbage, ",
                                " so we're closing the connection.")
                break
            else
                IPC_HANDLERS[msg_idx](channel, client_name, max_grid_byte_size)
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
function markovjunior_run_ipc(blocking::Bool
                              ;
                              max_grid_byte_size::Int=IPC_DEFAULT_MAX_GRID_BYTE_SIZE,
                              max_client_name_byte_size::Int=IPC_DEFAULT_MAX_CLIENT_NAME_BYTES,
                              allow_kill_message::Bool = false,
                              pipe_path::String = IPC_PIPE_PATH,
                              server_is_ready_callback = () -> println(stderr, "\tNow open to clients!"))
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
            max_grid_byte_size)
        )
    end
    if blocking
        server_loop()
        return nothing
    else
        Threads.@spawn server_loop()
        return server
    end
end

function markovjunior_run_ipc_main()::Cint
    arg_table = ArgParseSettings()
    @add_arg_table! arg_table begin
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
        "pipe_path"
            help = "An alternative name for the named-pipe that other processes talk to (valid syntax depends on your OS)"
            default = IPC_PIPE_PATH
    end

    try
        markovjunior_run_ipc(
            true,
            allow_kill_message = arg_table["allow-kill-msg"],
            max_grid_byte_size = arg_table["max-grid"],
            max_client_name_byte_size = arg_table["max-client-name"],
            pipe_path = arg_table["pipe_path"],
            server_is_ready_callback = () -> begin
                print(stderr, "\tNow open to clients! ")
                write(stdout, IPC_MAIN_START_CODE)
                println(stderr)
            end
        )
        println(stderr)
        write(stdout, IPC_MAIN_STOP_CODE)
        println(stderr)
        return 0
    catch e
        showerror(stderr, e, catch_backtrace())
        write(stdout, IPC_MAIN_STOP_CODE)
        return 1
    end
end
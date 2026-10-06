#include "jmj_connection.hpp"

#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>


using namespace jmj::ipc;

//Include named-pipe communication libraries.
//Sadly, while ASIO exists (and is already provided by Unreal engine),
//  it reportedly doesn't have the right abstractions for our use-case.
#if defined(_WIN32)

	#if JMJ_IS_UNREAL
		#include "Windows/AllowWindowsPlatformTypes.h"
		#include <Windows/WindowsHWrapper.h>
		#include "Windows/HideWindowsPlatformTypes.h"
	#else
        #ifndef WIN32_LEAN_AND_MEAN
            #define WIN32_LEAN_AND_MEAN
        #endif
        #ifndef NOMINMAX
		    #define NOMINMAX
        #endif
		#include <windows.h>
	#endif

	static_assert(sizeof(void*) == sizeof(HANDLE), "Can't cast between HANDLE and pimpl!");

    void* PipeConnection::GetNullPipeValue() { return static_cast<void*>(INVALID_HANDLE_VALUE); }

	void* PipeConnection::MakePipe()
	{
        //Windows uses UTF-16 for modern API calls; we need to convert the pipe name to that.
		std::wstring wPipeName;
		const auto& asciiPipeName = UsedPipeName();
		if (!asciiPipeName.empty())
		{
			const auto* asciiPipeNameBytes = reinterpret_cast<const char*>(asciiPipeName.data());
			int asciiPipeNameByteCount = static_cast<int>(asciiPipeName.size());
			DWORD conversionFlags = MB_ERR_INVALID_CHARS;
			
			int n = ::MultiByteToWideChar(CP_UTF8, conversionFlags, asciiPipeNameBytes,
										  asciiPipeNameByteCount, nullptr, 0);
			if (n <= 0)
			{
				InvokeErrorHandler(std::string("MultiByteToWideChar failed on pipe name *measurement*, with code ") + std::to_string(::GetLastError()));
				return static_cast<void*>(INVALID_HANDLE_VALUE);
			}

			wPipeName = std::wstring(static_cast<size_t>(n), L'\0'); //Note that brace-init would subtly pick a different constructor and cause a bug
			if (n != ::MultiByteToWideChar(CP_UTF8, conversionFlags, asciiPipeNameBytes,
							 		       asciiPipeNameByteCount, wPipeName.data(), n))
			{
				InvokeErrorHandler(std::string("MultiByteToWideChar failed on pipe name *conversion*, with code ") + std::to_string(::GetLastError()));
				return static_cast<void*>(INVALID_HANDLE_VALUE);
			}
		}
		
		while (true)
		{
			DWORD shareMode = 0; //Every call is blocking
			HANDLE hnd = ::CreateFileW(wPipeName.c_str(),
									   GENERIC_READ | GENERIC_WRITE, shareMode,
									   nullptr, OPEN_EXISTING, 0, nullptr);
			if (hnd != INVALID_HANDLE_VALUE)
				return static_cast<void*>(hnd);

			DWORD errCode = ::GetLastError();
			if (errCode == ERROR_FILE_NOT_FOUND)
			{
				InvokeErrorHandler(std::string("Pipe not found! You should have waited for the process to start up"));
				return static_cast<void*>(INVALID_HANDLE_VALUE);
			}
			else if (errCode == ERROR_PIPE_BUSY)
			{
				//Wait for the pipe to become available, then try again.
				::WaitNamedPipeW(wPipeName.c_str(), NMPWAIT_WAIT_FOREVER);
			}
			else
			{
				InvokeErrorHandler(std::string("Unexpected pipe connection failure on CreateFileW: ") + std::to_string(errCode));
				return static_cast<void*>(INVALID_HANDLE_VALUE);
			}
		}
	}
	void PipeConnection::ClosePipe(void*& pimpl)
	{
		auto& hnd = static_cast<HANDLE&>(pimpl);
		JMJ_ASSERT(hnd != INVALID_HANDLE_VALUE, "Tried to close a closed pipe!");
		
		::CloseHandle(hnd);
		hnd = INVALID_HANDLE_VALUE;
	}

	void PipeConnection::WriteToPipe(const std::byte* bytes, size_t count) const
	{
		HANDLE hnd = static_cast<HANDLE>(pipePimpl);
		CheckThread();

		size_t nWritten = 0;
		while (nWritten < count)
		{
			DWORD newWritten = 0;
            auto toWrite = static_cast<DWORD>(std::min(size_t{ MAXDWORD }, count - nWritten));
			if (!::WriteFile(hnd, bytes + nWritten, toWrite, &newWritten, nullptr))
			{
				InvokeErrorHandler(std::string("WriteFile to the JMJ IPC pipe failed, with code ") + std::to_string(::GetLastError()));
				return;
			}
			else if (newWritten < 1)
			{
				InvokeErrorHandler(std::string("Failed to write to the JMJ IPC pipe!"));
				return;
			}
			nWritten += newWritten;
		}
	}
	void PipeConnection::ReadFromPipe(std::byte* bytes, size_t count) const
	{
		HANDLE hnd = static_cast<HANDLE>(pipePimpl);
		CheckThread();

		size_t nRead = 0;
		while (nRead < count)
		{
			DWORD newNRead = 0;
            auto toRead = static_cast<DWORD>(std::min(size_t{ MAXDWORD }, count - nRead));
			if (!::ReadFile(hnd, bytes + nRead, toRead, &newNRead, nullptr))
			{
				InvokeErrorHandler(std::string("ReadFile to the JMJ IPC pipe failed, with code ") + std::to_string(::GetLastError()));
				return;
			}
			else if (newNRead < 1)
			{
				InvokeErrorHandler(std::string("ReadFile failed to read from the JMJ IPC pipe!"));
				return;
			}
			nRead += static_cast<size_t>(newNRead);
		}
	}


#elif defined(__linux__) || defined(__APPLE__)

	//IMPORTANT NOTE: this code was ported from the above Windows version by Claude;
	//  I use Windows so I can't test it.

	#include <sys/socket.h>
    #include <sys/un.h>
    #include <unistd.h>
    #include <fcntl.h>
    #include <cerrno>
    #include <cstring>
    #include <cstdint>

	static_assert(sizeof(void*) >= sizeof(int), "Can't cast between a socket fd and pimpl!");

    void* PipeConnection::GetNullPipeValue() { return reinterpret_cast<void*>(static_cast<std::intptr_t>(-1)); }

    namespace
    {
       int ToFd(void* pimpl) { return static_cast<int>(reinterpret_cast<std::intptr_t>(pimpl)); }
       void* FromFd(int fd) { return reinterpret_cast<void*>(static_cast<std::intptr_t>(fd)); }

       std::string ErrnoText(int e) { return std::to_string(e) + " (" + std::strerror(e) + ")"; }
    }

    void* PipeConnection::MakePipe()
    {
       //Unix socket paths are raw bytes, so ASCII needs no conversion -- only a length check.
       const auto& asciiPipeName = UsedPipeName();
       if (asciiPipeName.empty())
       {
          InvokeErrorHandler(std::string("The JMJ IPC pipe name is empty!"));
          return GetNullPipeValue();
       }

       sockaddr_un addr{ };
       addr.sun_family = AF_UNIX;
       //sun_path is only ~108 bytes (104 on macOS) and silently truncates if overrun.
       if (asciiPipeName.size() >= sizeof(addr.sun_path))
       {
          InvokeErrorHandler(std::string("The JMJ IPC pipe path is ") + std::to_string(asciiPipeName.size()) +
                               " bytes, but sun_path holds at most " + std::to_string(sizeof(addr.sun_path) - 1));
          return GetNullPipeValue();
       }
       std::memcpy(addr.sun_path, asciiPipeName.data(), asciiPipeName.size());

       //SOCK_CLOEXEC keeps this fd from leaking into any process we spawn later.
       #if defined(SOCK_CLOEXEC)
          int fd = ::socket(AF_UNIX, SOCK_STREAM | SOCK_CLOEXEC, 0);
       #else
          int fd = ::socket(AF_UNIX, SOCK_STREAM, 0);
          if (fd >= 0)
             ::fcntl(fd, F_SETFD, ::fcntl(fd, F_GETFD) | FD_CLOEXEC);
       #endif
       if (fd < 0)
       {
          InvokeErrorHandler(std::string("socket(AF_UNIX) failed, with code ") + ErrnoText(errno));
          return GetNullPipeValue();
       }

       //Without this, writing to a dead peer raises SIGPIPE and kills the process
       //   instead of returning an error the way Windows does.
       #if defined(SO_NOSIGPIPE)
          int on = 1;
          ::setsockopt(fd, SOL_SOCKET, SO_NOSIGPIPE, &on, sizeof(on));
       #endif

       while (true)
       {
          if (::connect(fd, reinterpret_cast<sockaddr*>(&addr), sizeof(addr)) == 0)
             return FromFd(fd);

          //Capture errno immediately; anything else we call can clobber it.
          int errCode = errno;
          if (errCode == EINTR)
             continue;

          ::close(fd);
          if (errCode == ENOENT)
             InvokeErrorHandler(std::string("Pipe not found! You should have waited for the process to start up"));
          else if (errCode == ECONNREFUSED)
             InvokeErrorHandler(std::string("The JMJ IPC socket file exists but nothing is accepting on it "
                                        "(stale file from a crashed run?)"));
          else
             InvokeErrorHandler(std::string("Unexpected pipe connection failure on connect(): ") + ErrnoText(errCode));
          return GetNullPipeValue();
       }
    }

    void PipeConnection::ClosePipe(void*& pimpl)
    {
       int fd = ToFd(pimpl);
       JMJ_ASSERT(fd != -1, "Tried to close a closed pipe!");

       //Never retry close() on EINTR; on Linux the descriptor is already gone.
       ::close(fd);
       pimpl = GetNullPipeValue();
    }

    void PipeConnection::WriteToPipe(const std::byte* bytes, size_t count) const
    {
       int fd = ToFd(pipePimpl);
       CheckThread();

       //MSG_NOSIGNAL is the Linux counterpart to SO_NOSIGPIPE above.
       #if defined(MSG_NOSIGNAL)
          constexpr int sendFlags = MSG_NOSIGNAL;
       #else
          constexpr int sendFlags = 0;
       #endif

       size_t written = 0;
       while (written < count)
       {
          ssize_t newWritten = ::send(fd, bytes + written, count - written, sendFlags);
          if (newWritten < 0)
          {
             int errCode = errno;
             if (errCode == EINTR)
                continue;
             InvokeErrorHandler(std::string("send() to the JMJ IPC pipe failed, with code ") + ErrnoText(errCode));
             return;
          }
          else if (newWritten < 1)
          {
             InvokeErrorHandler(std::string("Failed to write to the JMJ IPC pipe!"));
             return;
          }
          written += newWritten;
       }
    }

    void PipeConnection::ReadFromPipe(std::byte* bytes, size_t count) const
    {
       int fd = ToFd(pipePimpl);
       CheckThread();

       size_t nRead = 0;
       while (nRead < count)
       {
          ssize_t newNRead = ::recv(fd, bytes + nRead, count - nRead, 0);
          if (newNRead < 0)
          {
             int errCode = errno;
             if (errCode == EINTR)
                continue;
             InvokeErrorHandler(std::string("recv() from the JMJ IPC pipe failed, with code ") + ErrnoText(errCode));
             return;
          }
          else if (newNRead < 1)
          {
             //0 means a clean EOF: the peer closed. Treated as a failure, matching Windows.
             InvokeErrorHandler(std::string("recv() failed to read from the JMJ IPC pipe!"));
             return;
          }
          nRead += static_cast<size_t>(newNRead);
       }
    }

#else
	#error "MarkovJunior.jl runs as a separate process, and that is not currently supported on this platform"
#endif


PipeConnection::PipeConnection()
{
	owningThreadId = std::this_thread::get_id();
    pipePimpl = MakePipe();
    if (pipePimpl == GetNullPipeValue())
        InvokeErrorHandler("Failed to connect to the pipe!");

	//Generate a unique name for this client.
    static std::atomic<int> clientIndex = 0;
    auto thisClientIndex = clientIndex.fetch_add(1);
    clientLoggingName = ClientLoggingName();
    clientLoggingName += u8": ";
    clientLoggingName += std::u8string{ reinterpret_cast<const char8_t*>(std::to_string(thisClientIndex).c_str()) };

    //Perform the client handshake.
    WritePipeString(clientLoggingName);
}
PipeConnection::~PipeConnection()
{
	ClosePipe(pipePimpl);
}

const PipeConnection& PipeConnection::GetOrMakeThreadLocal()
{
	thread_local PipeConnection pipe{ };
    return pipe;
}
void PipeConnection::CheckThread() const
{
	JMJ_ASSERT(std::this_thread::get_id() == owningThreadId, "Called from the wrong thread!");
}
std::function<void(const std::string&)>& PipeConnection::ErrorHandler()
{
    static std::function<void(const std::string&)> handler = [](const std::string& msg) {
		#if JMJ_IS_UNREAL
    		FString msgU{ msg.c_str() };
    		UE_LOG(LogTemp, Fatal, TEXT("%s"), *msgU);
		#elif defined(__cpp_exceptions) || defined(__EXCEPTIONS) || defined(_CPPUNWIND) 
    		throw std::runtime_error{ msg.c_str() };
		#else
    		fprintf(stderr, "%s\n", msg.c_str());
    		std::abort();
		#endif
    };
    return handler;
}


std::variant<Handle_t, std::u8string> PipeConnection::ParseAlgorithm(std::string_view srcAscii) const
{
    thread_local std::u8string srcUtf8;
    srcUtf8.assign(reinterpret_cast<const char8_t*>(srcAscii.data()), srcAscii.size());
    return ParseAlgorithm(srcUtf8);
}
#if JMJ_IS_UNREAL
    std::variant<Handle_t, FString> PipeConnection::ParseAlgorithm(const FStringView srcUnreal) const
    {
        auto casted = StringCast<UTF8CHAR>(srcUnreal.GetData(), srcUnreal.Len());
        auto result = ParseAlgorithm(std::u8string_view{
            casted.Get(),
            static_cast<size_t>(casted.Length())
        });

        if (std::holds_alternative<Handle_t>(result))
            return std::get<Handle_t>(result);

        const auto& resultMsgUtf8 = std::get<std::u8string>(result);
        return FString::ConstructFromPtrSize(
            resultMsgUtf8.data(),
            static_cast<int>(resultMsgUtf8.size())
        );
    }
#endif
std::variant<Handle_t, std::u8string> PipeConnection::ParseAlgorithm(std::u8string_view srcUtf8) const
{
    CheckThread();

	WriteToPipe<uint32_t>(1);
	WritePipeString(srcUtf8);

    if (ReadPipeSuccessFlag())
    {
	    return ReadFromPipe<Handle_t>();
    }
    else
    {
    	std::u8string outErrorMsg;
	    ReadPipeString(outErrorMsg);
    	return outErrorMsg;
    }
}
bool PipeConnection::DestroyAlgorithm(Handle_t id) const
{
    CheckThread();

	WriteToPipe<uint32_t>(2);
    WriteToPipe(id);
    return ReadPipeSuccessFlag();
}

std::optional<Handle_t> PipeConnection::StartAlgorithmRun(Handle_t algo,
														  std::span<const int> resolution, std::span<const std::byte> seeds,
														  const TickSettings& tickSettings,
														  std::span<const Pixel_t> initialState) const
{
    CheckThread();

    size_t nElements = 1;
    for (int i : resolution)
        nElements *= i;
    if (!initialState.empty() && initialState.size() < nElements)
        return std::nullopt;

	WriteToPipe<uint32_t>(3);
    WriteToPipe(algo);

    //Write the desired resolution, and read if it was a valid value.
    WriteToPipe<uint32_t>(resolution.size());
    for (int s : resolution)
        WriteToPipe<uint32_t>(s);
    if (!ReadPipeSuccessFlag())
        return std::nullopt;

    //Write the initial state, if one was provided by the user.
    WriteToPipe<uint8_t>(initialState.empty() ? 0 : 1);
    if (!initialState.empty())
        WriteBytesToPipe(std::span{ initialState.data(), nElements });

    //Write the seeds.
    WriteToPipe<uint32_t>(seeds.size());
    WriteToPipe(seeds.data(), seeds.size());

    //Write tick settings.
    WriteToPipe<uint32_t>(tickSettings.MinPriorityLevel);
    WriteToPipe<uint8_t>(tickSettings.IsAnimationRender ? 1 : 0);

    //Get the result.
    if (ReadPipeSuccessFlag())
        return ReadFromPipe<Handle_t>();
    else
        return std::nullopt;
}
bool PipeConnection::CancelAlgorithmRun(Handle_t state) const
{
    CheckThread();

    WriteToPipe<uint32_t>(4);
    WriteToPipe<uint32_t>(state);
    return ReadPipeSuccessFlag();
}

TickResult PipeConnection::TickAlgorithm(Handle_t state,
                                         std::optional<TickAmount> amount,
                                         std::u8string& taggedEventNameBuffer) const
{
    CheckThread();

    WriteToPipe<uint32_t>(5);
    WriteToPipe(state);

    if (amount.has_value())
    {
        WriteToPipe<uint8_t>(0);
        WriteToPipe<uint32_t>(amount->PriorityLevel);
        WriteToPipe<uint32_t>(amount->Count);
    }
    else
    {
        WriteToPipe<uint8_t>(1);
    }

    if (!ReadPipeSuccessFlag())
        return TickResult::Error;
    if (ReadPipeSuccessFlag())
        return TickResult::AlgoCompleted;
    if (!ReadPipeSuccessFlag())
        return TickResult::Success;

    ReadPipeString(taggedEventNameBuffer);
    return TickResult::TaggedEvent;
}

bool PipeConnection::GetCurrentResolution(Handle_t state, std::vector<int>& inoutRes) const
{
    CheckThread();

    WriteToPipe<uint32_t>(11);
    WriteToPipe(state);
    if (!ReadPipeSuccessFlag())
        return false;

    auto nDims = ReadFromPipe<uint32_t>();
    for (decltype(nDims) i = 0; i < nDims; ++i)
        inoutRes.push_back(static_cast<int>(ReadFromPipe<uint32_t>()));
    return true;
}
bool PipeConnection::ReadCurrentGrid(Handle_t state, std::span<Pixel_t> outBuffer) const
{
    CheckThread();

    WriteToPipe<uint32_t>(6);
    WriteToPipe(state);
    if (!ReadPipeSuccessFlag())
        return false;

    auto nDims = ReadFromPipe<uint32_t>();
    size_t nElements = 1;
    for (decltype(nDims) i = 0; i < nDims; ++i)
        nElements *= ReadFromPipe<uint32_t>();
    if (outBuffer.size() < nElements)
    {
        //Drain the grid data from the pipe.
        //This is a waste of bandwidth, but it should be rare in practice.
        std::byte discardBuffer[4096];
        size_t remaining = nElements;
        while (remaining > 0)
        {
            size_t n = std::min(remaining, sizeof(discardBuffer));
            ReadFromPipe(discardBuffer, n);
            remaining -= n;
        }
        return false;
    }

    ReadBytesFromPipe(std::span{ outBuffer.data(), nElements });
    return true;
}
bool PipeConnection::WriteCurrentGrid(Handle_t state, std::span<const Pixel_t> srcBuffer) const
{
    CheckThread();

    WriteToPipe<uint32_t>(10);
    WriteToPipe(state);
    WriteToPipe<uint32_t>(srcBuffer.size() * sizeof(Pixel_t));

    if (!ReadPipeSuccessFlag())
        return false;
    WriteBytesToPipe(srcBuffer);

    auto successCode = ReadFromPipe<uint8_t>();
    JMJ_ASSERT(successCode == 1, "Expected code 1 but got %i", static_cast<int>(successCode));
    return true;
}
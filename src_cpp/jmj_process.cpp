#include "jmj_process.hpp"

#include <cstring>
#include <cstdio>
#include <algorithm>


using namespace jmj::ipc;

//Platform headers.
#if JMJ_IS_UNREAL
	#include "HAL/PlatformProcess.h"
#endif
#if defined(_WIN32)
	#if JMJ_IS_UNREAL
		#include "Windows/AllowWindowsPlatformTypes.h"
		#include <windows.h>
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
#elif defined(__linux__) || defined(__APPLE__)
	#include <sys/socket.h>
	#include <sys/un.h>
	#include <sys/wait.h>
	#include <unistd.h>
	#include <fcntl.h>
	#include <signal.h>
	#include <cerrno>
	#if defined(__linux__)
		#include <sys/prctl.h>
	#endif
#else
	#error "MarkovJunior.jl runs as a separate process, and that is not currently supported on this platform"
#endif


namespace
{
	//Adds one argument to a command line, quoted the way Windows' CommandLineToArgvW
	//   (and Unreal's own argument splitting) expects.
	template<typename String_t>
	void AppendQuotedArg(String_t& cmdLine, const String_t& arg)
	{
		using Char_t = typename String_t::value_type;
		if (!cmdLine.empty())
			cmdLine += Char_t(' ');

		bool needsQuotes = arg.empty() || std::any_of(arg.begin(), arg.end(), [](Char_t c) {
			return c == Char_t(' ') || c == Char_t('\t') || c == Char_t('\n') || c == Char_t('\v') || c == Char_t('"');
		});
		if (!needsQuotes)
		{
			cmdLine += arg;
			return;
		}

		//Backslashes are only special when they precede a quote.
		cmdLine += Char_t('"');
		size_t nBackslashes = 0;
		for (Char_t c : arg)
		{
			if (c == Char_t('\\'))
			{
				++nBackslashes;
				continue;
			}
			cmdLine.append((c == Char_t('"')) ? (nBackslashes * 2 + 1) : nBackslashes, Char_t('\\'));
			nBackslashes = 0;
			cmdLine += c;
		}
		cmdLine.append(nBackslashes * 2, Char_t('\\'));
		cmdLine += Char_t('"');
	}

	#if defined(_WIN32)
		std::wstring Utf8ToWide(std::string_view utf8)
		{
			if (utf8.empty())
				return { };
			int n = ::MultiByteToWideChar(CP_UTF8, 0, utf8.data(), static_cast<int>(utf8.size()), nullptr, 0);
			std::wstring wide(static_cast<size_t>(n), L'\0');
			::MultiByteToWideChar(CP_UTF8, 0, utf8.data(), static_cast<int>(utf8.size()), wide.data(), n);
			return wide;
		}

		//A job that kills its processes once its last handle closes,
		//   which the OS does for us if our own process dies.
		HANDLE MakeKillOnCloseJob(HANDLE hProcess)
		{
			HANDLE hJob = ::CreateJobObjectW(nullptr, nullptr);
			if (hJob == nullptr)
				return nullptr;

			JOBOBJECT_EXTENDED_LIMIT_INFORMATION limits{ };
			limits.BasicLimitInformation.LimitFlags = JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE;
			if (!::SetInformationJobObject(hJob, JobObjectExtendedLimitInformation, &limits, sizeof(limits)) ||
				!::AssignProcessToJobObject(hJob, hProcess))
			{
				::CloseHandle(hJob);
				return nullptr;
			}
			return hJob;
		}
	#endif

	//Whether some IPC server is already accepting clients on the given pipe.
	bool IsServerListening(const std::string& pipeName)
	{
		#if defined(_WIN32)
			//Doesn't consume a pipe instance, unlike actually connecting.
			std::wstring wideName = Utf8ToWide(pipeName);
			if (::WaitNamedPipeW(wideName.c_str(), 1))
				return true;
			return ::GetLastError() == ERROR_SEM_TIMEOUT; //Exists, but every instance is busy
		#else
			//A socket file can outlive a crashed server, so actually try connecting.
			//Note that the server will log this as a client failing its handshake.
			sockaddr_un addr{ };
			addr.sun_family = AF_UNIX;
			if (pipeName.empty() || pipeName.size() >= sizeof(addr.sun_path))
				return false;
			std::memcpy(addr.sun_path, pipeName.data(), pipeName.size());

			int fd = ::socket(AF_UNIX, SOCK_STREAM, 0);
			if (fd < 0)
				return false;
			bool connected = (::connect(fd, reinterpret_cast<sockaddr*>(&addr), sizeof(addr)) == 0);
			::close(fd);
			return connected;
		#endif
	}


	//The OS-specific half of the process manager.
	struct PlatformProcess
	{
		//Returns null and fills in the error message if the launch failed.
		static PlatformProcess* Launch(const ProcessSettings& settings, std::string& outError);
		//Closes our handles but doesn't kill the process (the kill-on-close job aside).
		~PlatformProcess();

		//Appends whatever output is available right now, without blocking.
		void ReadStdout(std::vector<char>& output) { ReadPipe(stdoutRead, output); }
		void ReadStderr(std::vector<char>& output) { ReadPipe(stderrRead, output); }

		bool IsRunning();
		//Asks the OS to kill the process; doesn't wait for it.
		void Terminate();

	#if JMJ_IS_UNREAL
		FProcHandle process;
		void* stdoutRead = nullptr;
		void* stderrRead = nullptr;
		#if PLATFORM_WINDOWS
			HANDLE killJob = nullptr;
		#endif
		TArray<uint8> readBuffer;

		void ReadPipe(void* pipe, std::vector<char>& output)
		{
			readBuffer.Reset();
			if (pipe != nullptr && FPlatformProcess::ReadPipeToArray(pipe, readBuffer))
				output.insert(output.end(), readBuffer.GetData(), readBuffer.GetData() + readBuffer.Num());
		}
	#elif defined(_WIN32)
		HANDLE process = nullptr;
		HANDLE stdoutRead = nullptr;
		HANDLE stderrRead = nullptr;
		HANDLE killJob = nullptr;

		static void ReadPipe(HANDLE pipe, std::vector<char>& output)
		{
			//Anonymous pipes can't be non-blocking, so only ask for what's already there.
			DWORD available = 0;
			if (!::PeekNamedPipe(pipe, nullptr, 0, nullptr, &available, nullptr) || available == 0)
				return;

			size_t start = output.size();
			output.resize(start + available);
			DWORD nRead = 0;
			if (!::ReadFile(pipe, output.data() + start, available, &nRead, nullptr))
				nRead = 0;
			output.resize(start + nRead);
		}
	#else
		pid_t pid = -1;
		int stdoutRead = -1;
		int stderrRead = -1;
		bool reaped = false;

		static void ReadPipe(int fd, std::vector<char>& output)
		{
			char buffer[4096];
			while (true)
			{
				ssize_t n = ::read(fd, buffer, sizeof(buffer));
				if (n > 0)
					output.insert(output.end(), buffer, buffer + n);
				else if (n < 0 && errno == EINTR)
					continue;
				else
					return; //0 is end-of-stream, EAGAIN is "nothing yet"
			}
		}
	#endif
	};


#if JMJ_IS_UNREAL

	PlatformProcess* PlatformProcess::Launch(const ProcessSettings& settings, std::string& outError)
	{
		auto* proc = new PlatformProcess();

		//Unreal takes the arguments as one string.
		std::basic_string<TCHAR> cmdLine;
		for (const FString& arg : settings.Arguments)
			AppendQuotedArg(cmdLine, std::basic_string<TCHAR>{ *arg, static_cast<size_t>(arg.Len()) });

		void *stdoutWrite = nullptr,
			 *stderrWrite = nullptr;
		if (!FPlatformProcess::CreatePipe(proc->stdoutRead, stdoutWrite) ||
			!FPlatformProcess::CreatePipe(proc->stderrRead, stderrWrite))
		{
			outError = "Couldn't create the stdout/stderr pipes";
			FPlatformProcess::ClosePipe(proc->stdoutRead, stdoutWrite);
			proc->stdoutRead = nullptr;
			delete proc;
			return nullptr;
		}

		proc->process = FPlatformProcess::CreateProc(*settings.ExecutablePath, cmdLine.c_str(),
													 false, settings.Hidden, settings.Hidden,
													 nullptr, 0, nullptr,
													 stdoutWrite, nullptr, stderrWrite);
		//The child has its own copies of the write ends now.
		FPlatformProcess::ClosePipe(nullptr, stdoutWrite);
		FPlatformProcess::ClosePipe(nullptr, stderrWrite);
		if (!proc->process.IsValid())
		{
			auto exePathUtf8 = StringCast<UTF8CHAR>(*settings.ExecutablePath);
			outError = "CreateProc failed for " + std::string(reinterpret_cast<const char*>(exePathUtf8.Get()),
															  static_cast<size_t>(exePathUtf8.Length()));
			delete proc;
			return nullptr;
		}

		#if PLATFORM_WINDOWS
			//Unreal can't launch suspended, so the process runs briefly before joining the job.
			proc->killJob = MakeKillOnCloseJob(proc->process.Get());
		#endif
		return proc;
	}
	PlatformProcess::~PlatformProcess()
	{
		FPlatformProcess::ClosePipe(stdoutRead, nullptr);
		FPlatformProcess::ClosePipe(stderrRead, nullptr);
		if (process.IsValid())
			FPlatformProcess::CloseProc(process);
		#if PLATFORM_WINDOWS
			if (killJob != nullptr)
				::CloseHandle(killJob);
		#endif
	}
	bool PlatformProcess::IsRunning() { return FPlatformProcess::IsProcRunning(process); }
	void PlatformProcess::Terminate() { FPlatformProcess::TerminateProc(process, true); }

#elif defined(_WIN32)

	PlatformProcess* PlatformProcess::Launch(const ProcessSettings& settings, std::string& outError)
	{
		auto* proc = new PlatformProcess();
		auto fail = [&](const char* what) {
			outError = std::string(what) + " failed, with code " + std::to_string(::GetLastError());
			delete proc;
			return nullptr;
		};

		std::wstring exePath = Utf8ToWide(settings.ExecutablePath);
		std::wstring cmdLine;
		AppendQuotedArg(cmdLine, exePath);
		for (const auto& arg : settings.Arguments)
			AppendQuotedArg(cmdLine, Utf8ToWide(arg));

		//The child's ends of the pipes must be inheritable, and ours must not be.
		SECURITY_ATTRIBUTES inheritable{ sizeof(SECURITY_ATTRIBUTES), nullptr, TRUE };
		HANDLE stdoutWrite = nullptr,
			   stderrWrite = nullptr;
		if (!::CreatePipe(&proc->stdoutRead, &stdoutWrite, &inheritable, 0))
			return fail("CreatePipe(stdout)");
		if (!::CreatePipe(&proc->stderrRead, &stderrWrite, &inheritable, 0))
		{
			::CloseHandle(stdoutWrite);
			return fail("CreatePipe(stderr)");
		}
		::SetHandleInformation(proc->stdoutRead, HANDLE_FLAG_INHERIT, 0);
		::SetHandleInformation(proc->stderrRead, HANDLE_FLAG_INHERIT, 0);
		HANDLE nul = ::CreateFileW(L"NUL", GENERIC_READ, FILE_SHARE_READ | FILE_SHARE_WRITE, &inheritable,
								   OPEN_EXISTING, 0, nullptr);

		STARTUPINFOW startup{ };
		startup.cb = sizeof(startup);
		startup.dwFlags = STARTF_USESTDHANDLES;
		startup.hStdInput = nul;
		startup.hStdOutput = stdoutWrite;
		startup.hStdError = stderrWrite;

		//Start suspended, so it's in the kill-on-close job before it can do anything.
		DWORD flags = CREATE_SUSPENDED | (settings.Hidden ? CREATE_NO_WINDOW : 0);
		PROCESS_INFORMATION info{ };
		BOOL launched = ::CreateProcessW(exePath.c_str(), cmdLine.data(), nullptr, nullptr, TRUE,
										 flags, nullptr, nullptr, &startup, &info);
		DWORD launchError = ::GetLastError();
		::CloseHandle(stdoutWrite);
		::CloseHandle(stderrWrite);
		if (nul != INVALID_HANDLE_VALUE)
			::CloseHandle(nul);
		if (!launched)
		{
			::SetLastError(launchError);
			return fail("CreateProcessW");
		}

		proc->process = info.hProcess;
		proc->killJob = MakeKillOnCloseJob(info.hProcess);
		::ResumeThread(info.hThread);
		::CloseHandle(info.hThread);
		return proc;
	}
	PlatformProcess::~PlatformProcess()
	{
		for (HANDLE h : { stdoutRead, stderrRead, process, killJob })
			if (h != nullptr)
				::CloseHandle(h);
	}
	bool PlatformProcess::IsRunning() { return ::WaitForSingleObject(process, 0) == WAIT_TIMEOUT; }
	void PlatformProcess::Terminate() { ::TerminateProcess(process, 1); }

#else

	PlatformProcess* PlatformProcess::Launch(const ProcessSettings& settings, std::string& outError)
	{
		//Prepare everything before forking; the child may only make async-signal-safe calls.
		std::vector<std::string> argStorage{ settings.ExecutablePath };
		argStorage.insert(argStorage.end(), settings.Arguments.begin(), settings.Arguments.end());
		std::vector<char*> argv;
		for (auto& arg : argStorage)
			argv.push_back(arg.data());
		argv.push_back(nullptr);

		int stdoutPipe[2], stderrPipe[2];
		if (::pipe(stdoutPipe) != 0)
		{
			outError = "pipe() failed, with errno " + std::to_string(errno);
			return nullptr;
		}
		if (::pipe(stderrPipe) != 0)
		{
			outError = "pipe() failed, with errno " + std::to_string(errno);
			::close(stdoutPipe[0]);
			::close(stdoutPipe[1]);
			return nullptr;
		}
		//Keep these out of any other process we spawn; dup2() below clears it for the child's copies.
		for (int fd : { stdoutPipe[0], stdoutPipe[1], stderrPipe[0], stderrPipe[1] })
			::fcntl(fd, F_SETFD, FD_CLOEXEC);

		pid_t parentPid = ::getpid();
		pid_t pid = ::fork();
		if (pid == 0)
		{
			::dup2(stdoutPipe[1], STDOUT_FILENO);
			::dup2(stderrPipe[1], STDERR_FILENO);
			int devNull = ::open("/dev/null", O_RDONLY);
			if (devNull >= 0)
				::dup2(devNull, STDIN_FILENO);

			#if defined(__linux__)
				//Die with the parent. Linux ties this to the parent *thread*,
				//   which is why ProcessManager launches us from its monitor thread.
				::prctl(PR_SET_PDEATHSIG, SIGKILL);
				if (::getppid() != parentPid)
					::_exit(127);
			#endif

			::execv(argv[0], argv.data());
			::_exit(127);
		}

		::close(stdoutPipe[1]);
		::close(stderrPipe[1]);
		if (pid < 0)
		{
			outError = "fork() failed, with errno " + std::to_string(errno);
			::close(stdoutPipe[0]);
			::close(stderrPipe[0]);
			return nullptr;
		}

		auto* proc = new PlatformProcess();
		proc->pid = pid;
		proc->stdoutRead = stdoutPipe[0];
		proc->stderrRead = stderrPipe[0];
		for (int fd : { proc->stdoutRead, proc->stderrRead })
			::fcntl(fd, F_SETFL, ::fcntl(fd, F_GETFL) | O_NONBLOCK);
		return proc;
	}
	PlatformProcess::~PlatformProcess()
	{
		::close(stdoutRead);
		::close(stderrRead);
		if (!reaped)
			::waitpid(pid, nullptr, WNOHANG);
	}
	bool PlatformProcess::IsRunning()
	{
		if (reaped)
			return false;
		int status = 0;
		pid_t result = ::waitpid(pid, &status, WNOHANG);
		if (result == 0)
			return true;
		reaped = true;
		return false;
	}
	void PlatformProcess::Terminate() { if (!reaped) ::kill(pid, SIGKILL); }

#endif
}


ProcessManager::ProcessManager(ProcessSettings _settings)
	: settings(std::move(_settings))
{
	if (!settings.StderrLineHandler)
	{
		settings.StderrLineHandler = [](ProcessStringView_t line) {
			#if JMJ_IS_UNREAL
				UE_LOG(LogTemp, Log, TEXT("<JMJ IPC> %.*s"), line.Len(), line.GetData());
			#else
				std::fprintf(stderr, "<JMJ IPC> %.*s\n", static_cast<int>(line.size()), line.data());
			#endif
		};
	}
	auto log = [&](const std::string& msg) { EmitLine(msg); };

	if (settings.ReuseExistingServer && IsServerListening(PipeConnection::UsedPipeName()))
	{
		log("[client] Found an existing IPC server, so we won't launch our own");
		isManagingProcess = false;
		SetState(ProcessState::Ready);
		return;
	}

	//The monitor thread launches the process itself.
	isManagingProcess = true;
	SetState(ProcessState::Booting);
	monitorThread = std::thread([this]() { MonitorLoop(); });
}
ProcessManager::~ProcessManager()
{
	stopMonitor = true;
	if (monitorThread.joinable())
		monitorThread.join();

	if (auto* proc = static_cast<PlatformProcess*>(platformPimpl.load()))
	{
		if (proc->IsRunning())
			proc->Terminate();
		delete proc;
	}
}

void ProcessManager::SetState(ProcessState newState)
{
	{
		std::lock_guard lock{ stateMutex };
		state = newState;
	}
	stateChanged.notify_all();
}
bool ProcessManager::WaitUntilReady(std::optional<std::chrono::milliseconds> timeout) const
{
	std::unique_lock lock{ stateMutex };
	auto isSettled = [&]() { return state.load() != ProcessState::Booting; };
	if (timeout.has_value())
		stateChanged.wait_for(lock, *timeout, isSettled);
	else
		stateChanged.wait(lock, isSettled);
	return state.load() == ProcessState::Ready;
}
void ProcessManager::Kill()
{
	//The monitor thread notices the death and updates the state.
	//If the process is still being launched, the monitor thread kills it right after.
	killRequested = true;
	if (auto* proc = static_cast<PlatformProcess*>(platformPimpl.load()))
		proc->Terminate();
}

void ProcessManager::EmitLine(std::string_view utf8) const
{
	#if JMJ_IS_UNREAL
		FString line = FString::ConstructFromPtrSize(reinterpret_cast<const UTF8CHAR*>(utf8.data()),
													 static_cast<int32>(utf8.size()));
		settings.StderrLineHandler(line);
	#else
		settings.StderrLineHandler(utf8);
	#endif
}

void ProcessManager::MonitorLoop()
{
	//Launch from this thread, which lives exactly as long as we want the process to.
	std::string launchError;
	auto* proc = PlatformProcess::Launch(settings, launchError);
	if (proc == nullptr)
	{
		EmitLine("[client] Failed to launch the IPC process: " + launchError);
		SetState(ProcessState::Dead);
		return;
	}
	platformPimpl = proc;
	if (killRequested.load())
		proc->Terminate();

	std::vector<char> stdoutBytes, stderrBytes;

	//stdout only carries the 4-byte start/stop codes.
	auto handleStdoutCodes = [&]()
	{
		while (stdoutBytes.size() >= sizeof(uint32_t))
		{
			uint32_t code;
			std::memcpy(&code, stdoutBytes.data(), sizeof(code));
			stdoutBytes.erase(stdoutBytes.begin(), stdoutBytes.begin() + sizeof(code));

			if (code == StdoutStartCode && state.load() == ProcessState::Booting)
				SetState(ProcessState::Ready);
			else if (code == StdoutStopCode)
				SetState(ProcessState::Dying);
			else
				EmitLine("[client] Unexpected code on the IPC process's stdout: " + std::to_string(code));
		}
	};
	//Forward complete lines of stderr (or everything, once the process is gone).
	auto handleStderrLines = [&](bool includePartialLine)
	{
		size_t lineStart = 0;
		for (size_t i = 0; i < stderrBytes.size(); ++i)
		{
			if (stderrBytes[i] != '\n')
				continue;
			size_t lineEnd = (i > lineStart && stderrBytes[i - 1] == '\r') ? (i - 1) : i;
			EmitLine({ stderrBytes.data() + lineStart, lineEnd - lineStart });
			lineStart = i + 1;
		}
		if (includePartialLine && lineStart < stderrBytes.size())
		{
			EmitLine({ stderrBytes.data() + lineStart, stderrBytes.size() - lineStart });
			lineStart = stderrBytes.size();
		}
		stderrBytes.erase(stderrBytes.begin(), stderrBytes.begin() + lineStart);
	};

	while (!stopMonitor.load())
	{
		//Check for death before reading, so the final read below can't miss output.
		bool isRunning = proc->IsRunning();

		proc->ReadStdout(stdoutBytes);
		proc->ReadStderr(stderrBytes);
		handleStdoutCodes();
		handleStderrLines(!isRunning);

		if (!isRunning)
		{
			SetState(ProcessState::Dead);
			return;
		}

		//The pipes are small, so poll often enough that a chatty process never blocks on stderr.
		std::this_thread::sleep_for(std::chrono::milliseconds(20));
	}
}

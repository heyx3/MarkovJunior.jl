#pragma once

#include <string>
#include <string_view>
#include <vector>
#include <functional>
#include <atomic>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <chrono>
#include <optional>

//For JMJ_IS_UNREAL, JMJ_MODULE_EXPORT, JMJ_ASSERT, and the pipe name.
#include "jmj_connection.hpp"


namespace jmj
{
	namespace ipc
	{
		//Conventional string types for process management:
		//  Unreal's in Unreal projects, otherwise the standard ones (assumed to hold UTF-8).
		#if JMJ_IS_UNREAL
			using ProcessString_t = FString;
			using ProcessStringView_t = FStringView;
			using ProcessStringList_t = TArray<FString>;
		#else
			using ProcessString_t = std::string;
			using ProcessStringView_t = std::string_view;
			using ProcessStringList_t = std::vector<std::string>;
		#endif

		//The lifecycle of the IPC process, as seen from our side.
		enum class ProcessState
		{
			//Launched, but not accepting clients yet.
			Booting,
			//Accepting clients.
			Ready,
			//No longer accepting new clients; it exits once the existing ones disconnect.
			Dying,
			//Not running (it exited, failed to launch, or was killed).
			Dead
		};

		struct ProcessSettings
		{
			//Path to the IPC executable (e.g. ".../bin/JMarkovJunior_IPC.exe").
			ProcessString_t ExecutablePath;
			//Command-line arguments, passed through as-is (quoting is handled for you).
			//Remember that the pipe name is the positional argument, if you changed 'PipeConnection::UsedPipeName()'.
			ProcessStringList_t Arguments;

			//If true, the process gets no console window (only matters on Windows).
			bool Hidden = true;

            //If true, checks for an existing IPC pipe before creating one.
            //This allows you to run it from a Julia script/REPL.
			bool ReuseExistingServer = true;

			//Receives each line of the IPC process's stderr (without the newline).
			//Called from a background thread!
			//Defaults to Unreal logging in Unreal projects, otherwise pipes it through to stderr.
			std::function<void(ProcessStringView_t)> StderrLineHandler;
		};

		//Launches and watches the IPC process on a background thread.
		//The process is killed when this object is destroyed,
		//   and on Windows/Linux the OS also kills it if our own process dies first.
		//Destroy this only after you're done with every 'PipeConnection',
		//   otherwise their next call hits the fatal-error handler.
		class JMJ_MODULE_EXPORT ProcessManager
		{
		public:

			explicit ProcessManager(ProcessSettings settings);
			~ProcessManager();

			ProcessManager(const ProcessManager&) = delete;
			ProcessManager& operator=(const ProcessManager&) = delete;

			ProcessState GetState() const { return state.load(); }
			//False if we're connected to an IPC process we didn't start (usually done for debugging).
			bool IsManagingProcess() const { return isManagingProcess; }

			//Blocks until the process is ready for clients, it dies, or the timeout passes.
			//Returns whether it's ready.
			bool WaitUntilReady(std::optional<std::chrono::milliseconds> timeout = std::nullopt) const;

			//Forcibly kills the process IF we own it.
			void Kill();

		private:

			ProcessSettings settings;
			bool isManagingProcess = false;

			std::atomic<ProcessState> state{ ProcessState::Dead };
			mutable std::mutex stateMutex;
			mutable std::condition_variable stateChanged;
			void SetState(ProcessState newState);

			//The process is launched *from* the monitor thread,
			//   because Linux ties the kill-on-parent-death signal to the launching thread.
			std::atomic<void*> platformPimpl{ nullptr };
			std::atomic<bool> stopMonitor{ false },
							  killRequested{ false };
			std::thread monitorThread;
			void MonitorLoop();

			//Sends UTF-8 text to the stderr handler.
			void EmitLine(std::string_view utf8) const;
		};
	}
}

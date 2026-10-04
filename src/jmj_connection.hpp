#pragma once

#include <span>
#include <variant>
#include <string>
#include <cstdint>
#include <thread>
#include <functional>
#include <optional>
#include <cassert>
#include <vector>

#include "jmj_consts.hpp"
#include "jmj_ipc.hpp"

//A lot of defaults (like the assert macro) will change from normal C++ to Unreal C++.
#if defined(UE_BUILD_DEBUG) || defined(UE_BUILD_DEVELOPMENT) || defined(UE_BUILD_TEST) || defined(UE_BUILD_SHIPPING)
	#define JMJ_IS_UNREAL 1
    #include "CoreMinimal.h"
#else
	#define JMJ_IS_UNREAL 0
#endif

#ifndef JMJ_ASSERT
	#if JMJ_IS_UNREAL
		#define JMJ_ASSERT(cond, msg, ...) checkf(cond, TEXT(msg) __VA_OPT__(,) __VA_ARGS__)
	#else
		#define JMJ_ASSERT(cond, msg, ...) assert(cond && msg)
	#endif
#endif

#ifndef JMJ_MODULE_EXPORT
    #define JMJ_MODULE_EXPORT
#endif


//MarkovJunior.jl
namespace jmj
{
	//The connection to the other process running MarkovJunior.jl ("IPC")
	namespace ipc
	{
		//The unsigned type representing an algorithm/state handle in the IPC.
		using Handle_t = uint32_t;
		//The unsigned type representing a grid color value.
		using Pixel_t = uint8_t;

	    //The initial settings when starting an algorithm run.
	    struct TickSettings
	    {
	        //Controls how granular the algorithm's ticks are;
	        //  on worker threads you should normally use high-priority to only count large jumps in progress.
	        //If generating animations, use a low value like 1 or 2!
	        int MinPriorityLevel = jmj::TickPriorityHigh;

	        //Disables certain optimizations that would change the order of events,
	        //   so that the state in between every tick looks correct when animated.
	        bool IsAnimationRender = false;
	    };
		//Possible outcomes from advancing an algorithm state.
		enum class TickResult
		{
		    //The tick completed without incident (no tagged events, no end of algorithm).
			Success,
		    //The tick did not run because your parameters were bad.
		    Error,
		    //A named event happened.
			TaggedEvent,
		    //The algorithm has finished.
		    AlgoCompleted
		};
		struct TickAmount
		{
			int PriorityLevel = jmj::TickPriorityHigh;
			int Count = 1;
		};

		//One instance of a connection to the IPC's named pipe.
		//Manages itself as a per-thread singleton.
		struct JMJ_MODULE_EXPORT PipeConnection
		{
			//The name of these connection clients to the IPC, only for logging purposes.
			//Each individual client appends a unique index to this name.
			static std::u8string& ClientLoggingName() { static std::u8string name{ u8"C++ Client" }; return name; }
			//The name of the pipe to connect to; defaults to the IPC's own default pipe name.
			//Users can change this *before* the IPC starts and the first connection is made!
			static std::string& UsedPipeName() { static std::string name{ jmj::ipc::NamedPipe }; return name; }
			//The callback to invoke when there is a fatal IPC error (e.g. the IPC process died).
			//Default behavior uses Unreal fatal-error logging in Unreal projects,
			//  or else an exception in exception-enabled code,
			//  or else prints to stderr and then exits.
			static std::function<void(const std::string&)>& ErrorHandler();

			//Gets (or sets up) an IPC connection for the calling thread.
			//All interaction with it must *stay* on this thread!
			static const PipeConnection& GetOrMakeThreadLocal();
			// (NOTE: must stay in the cpp or else different instances will exist across DLL boundaries)

			//Gets the unique name of this client in the IPC process's log.
			const auto& GetClientLoggingName() const { return clientLoggingName; }

			//Returns the handle of the new algorithm, or else an error message about your source code.
			std::variant<Handle_t, std::u8string> ParseAlgorithm(std::string_view srcAscii) const;
			//Returns the handle of the new algorithm, or else an error message about your source code.
			std::variant<Handle_t, std::u8string> ParseAlgorithm(std::u8string_view srcUtf8) const;
            #if JMJ_IS_UNREAL
			    //Returns the handle of the new algorithm, or else an error message about your source code.
                std::variant<Handle_t, FString> ParseAlgorithm(const FStringView srcUnreal) const;
            #endif
			//Does not affect ongoing algorithm runs.
			//Returns whether it was successful (i.e. if the handle was a real algorithm).
			bool DestroyAlgorithm(Handle_t) const;

			//Returns the handle of the new running algorithm state, or null if your parameters were bad.
			//You can customize the starting state of the algorithm (providing an array with X axis innermost),
			//  otherwise it'll start with its default fill color.
			std::optional<Handle_t> StartAlgorithmRun(Handle_t algo,
                                                      std::span<const int> resolution, std::span<const std::byte> seeds,
												      const TickSettings& tickSettings = { },
										  		      std::span<const Pixel_t> initialState = { }) const;
			//Returns whether it was successful (i.e. if the handle was a real algorithm state).
			bool CancelAlgorithmRun(Handle_t state) const;
			//Advances the algorithm by the given number of ticks, or until a custom named event (`@event`) happens.
			//Of course it will also stop if the algorithm completes.
			//
			//Returns what happened (and if it was a tagged event, writes its name into your provided string).
			//
			//Note that the tick priority you provide here is capped below
            //   by the priority you gave on algorithm start.
			TickResult TickAlgorithm(Handle_t state, std::optional<TickAmount> amount,
									 std::u8string& taggedEventNameBuffer) const;

			//Queries the given state's grid resolution, and appends it into the given array.
			//Returns whether this was successful (i.e. the ID was valid).
			bool GetCurrentResolution(Handle_t state, std::vector<int>& inoutRes) const;
			//Reads the current grid state into the given buffer.
			//Returns whether it succeeded (ID is good and buffer is big enough).
			bool ReadCurrentGrid(Handle_t state, std::span<Pixel_t> outBuffer) const;
			//Writes a new grid into the given algorithm state (only allowed at particular points, be careful!)
			//Returns whether it succeeded (ID is good and source buffer is correct size).
			bool WriteCurrentGrid(Handle_t state, std::span<const Pixel_t> srcBuffer) const;


		private:

			PipeConnection();
			~PipeConnection();

			//As a singleton, it doesn't need to go anywhere.
			PipeConnection(const PipeConnection&) = delete;
			PipeConnection(PipeConnection&&) = delete;

			std::u8string clientLoggingName;

			void* pipePimpl;
            static void* GetNullPipeValue();
			static void* MakePipe();
			static void ClosePipe(void*&);
			void WriteToPipe(const std::byte* bytes, size_t count) const;
			void ReadFromPipe(std::byte* bytes, size_t count) const;

			std::thread::id owningThreadId;
			void CheckThread() const;

            static void InvokeErrorHandler(const std::string& msg)
            {
                ErrorHandler()(msg);
                std::abort();
            }


			template<typename T>
			void WriteToPipe(const T& bitsData) const { WriteToPipe(reinterpret_cast<const std::byte*>(&bitsData), sizeof(T)); }
			template<typename T>
			void WriteBytesToPipe(std::span<const T> elements) const { WriteToPipe(reinterpret_cast<const std::byte*>(elements.data()), elements.size() * sizeof(T)); }

			template<typename T>
			void ReadFromPipe(T& output) const { ReadFromPipe(reinterpret_cast<std::byte*>(&output), sizeof(T)); }
			template<typename T>
			T ReadFromPipe() const { T output; ReadFromPipe(output); return output; }
			template<typename T>
			void ReadBytesFromPipe(std::span<T> output) const
			{
			    ReadFromPipe(reinterpret_cast<std::byte*>(output.data()),
			                 static_cast<size_t>(output.size()) * sizeof(T));
			}

			bool ReadPipeSuccessFlag() const { return ReadFromPipe<uint8_t>() == 1; }
			void ReadPipeString(std::u8string& outputBuffer) const
			{
				outputBuffer.resize(ReadFromPipe<uint32_t>());
				ReadFromPipe(reinterpret_cast<std::byte*>(outputBuffer.data()),
							 outputBuffer.size());

                JMJ_ASSERT(!outputBuffer.empty() && outputBuffer.back() == u8'\0',
                           "String should come with its null-terminator!");
                outputBuffer.pop_back();
			}
			void WritePipeString(std::span<const char8_t> utf8) const
			{
				WriteToPipe<uint32_t>(utf8.size());
				WriteToPipe(reinterpret_cast<const std::byte*>(utf8.data()),
                            utf8.size());
			}
		};
	}
}
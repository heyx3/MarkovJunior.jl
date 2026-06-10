#ifndef JMJ_LIB_INTERFACE_H
#define JMJ_LIB_INTERFACE_H

#ifdef __cplusplus
    extern "C" {
#endif


//Users of this header should deine JMJ_API as appropriate.
#ifndef JMJ_API
    #define JMJ_API
#endif

typedef unsigned int JmjID;

//Julia's standard init and shutdown functions; we copied them to here so we can handle their linkage.
JMJ_API void init_julia(int argc, char *argv[]);
JMJ_API void shutdown_julia(int retcode);

//Any time an error would be printed to stderr, you can instead drop it by setting this flag to true.
JMJ_API void jmj_set_suppress_stderr(int b);
JMJ_API int jmj_get_suppress_stderr();

//Parses the given UTF-8 text as a @markovjunior algorithm.
//If it succeeded, returns its ID (a positive sequential integer).
//Otherwise, writes into your error string buffer (also UTF-8) and returns 0.
//
//If you don't provide an error string buffer, the error is printed to stderr instead.
//
//When you're done using the parsed instance, you should call jmj_close_algo().
JMJ_API JmjID jmj_algo_parse(const char* algoStr, char* outErrStr, int errStrBufSize);
//Destroys an algorithm that you created with 'jmj_algo_parse()'.
JMJ_API void jmj_algo_close(JmjID algoID);

//Starts running an instance of the given algorithm,
//    on a grid of the given size,
//    and optionally using the given seed (or a random one if no seed is provided).
//
//Returns the ID of the new algorithm state (a positive sequential integer),
//    or 0 if your arguments were invalid.
JMJ_API JmjID jmj_start(JmjID algoID,
                        int nDims, const int* size,
                        int nSeedBytes, const unsigned char* seedBytes);
//Cleans up the given algorithm state, whether it was done running or not.
//This also invalidates the grid memory ('jmj_grid()').
JMJ_API void jmj_destroy(JmjID stateID);


//Runs some number of ticks of the given algorithm instance,
//    and returns whether it finished.
JMJ_API int jmj_step(JmjID algoID, JmjID stateID, int count);
//Runs the algorithm until it finishes.
//Beware of infinite loops!
JMJ_API void jmj_finish(JmjID algoID, JmjID stateID);
//Checks whether an algorithm instance has finished running.
JMJ_API int jmj_is_finished(JmjID algoID, JmjID stateI);

//Retrieves the grid that the given instance is operating on.
//This memory is only safe until the next tick, as some algorithm ops will reallocate the grid!
JMJ_API const unsigned char* jmj_grid(JmjID stateID);


#ifdef __cplusplus
    }
#endif

#endif
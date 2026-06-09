typedef unsigned int JmjID;

//Any time an error would be printed to stderr, you can instead drop it by setting this flag to true.
void jmj_set_suppress_stderr(int b);
int jmj_get_suppress_stderr();

//Parses the given text as a @markovjunior algorithm.
//If it succeeded, returns its ID (a positive sequential integer).
//Otherwise, writes into your error string and returns 0.
//
//If you don't provide an error string buffer, the error is printed to stderr instead.
//
//When you're done using the parsed instance, you should call jmj_close_algo().
JmjID jmj_algo_parse(const char* algoStr, char* outErrStr, int errStrBufSize);
//Destroys an algorithm that you created with 'jmj_algo_parse()'.
void jmj_algo_close(JmjID algoID);

//Starts running an instance of the given algorithm,
//    on a grid of the given size,
//    and optionally using the given seed (or a random one if no seed is provided).
//
//Returns the ID of the new algorithm state (a positive sequential integer),
//    or 0 if your arguments were invalid.
JmjID jmj_start(JmjID algoID,
                int nDims, const int* size,
                int nSeedBytes, const unsigned char* seedBytes);
//Checks whether an algorithm instance has finished running.
int jmj_is_finished(JmjID algoID, JmjID stateI);

//Runs some number of ticks of the given algorithm instance,
//    and returns whether it finished.
int jmj_step(JmjID algoID, JmjID stateID, int count);
//Runs the algorithm until it finishes.
//Beware of infinite loops!
void jmj_finish(JmjID algoID, JmjID stateID);

//Retrieves the grid that the given instance is operating on.
//This memory is only safe until the next tick, as some algorithm ops will reallocate the grid!
const unsigned char* jmj_grid(JmjID stateID);
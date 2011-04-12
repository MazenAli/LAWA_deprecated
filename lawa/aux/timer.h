#ifndef LAWA_ADAPTIVE_AUX_TIMER_H
#define LAWA_ADAPTIVE_AUX_TIMER_H 1

#include <ctime>

struct Timer
{
	void
	start()
	{
	    _start = clock();
	}

    void
    stop()
    {
        _stop = clock();
    }

    double
    elapsed()
    {
        return double(_stop - _start) / CLOCKS_PER_SEC;
    }

    clock_t _start;
    clock_t _stop;
};



#endif //LAWA_ADAPTIVE_AUX_TIMER_H 1

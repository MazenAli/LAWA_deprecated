#ifndef LAWA_AUX_TIMER_H
#define LAWA_AUX_TIMER_H 1

#include <ctime>

struct Timer
{
	void
    start();

    void
    stop();

    double
    elapsed();

    clock_t _start;
    clock_t _stop;
};

#endif //LAWA_AUX_TIMER_H

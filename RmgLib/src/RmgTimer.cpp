/*
 *
 * Copyright (c) 2014, Emil Briggs
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.

 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
*/

#if !(defined(_WIN32) || defined(_WIN64))
    #include <sys/time.h>
#endif
#include <stdio.h>
#include <string>
#include <iostream>
#include <functional>
#include <string>
#include <iomanip>
#include <ctime>
#include <ratio>
#include <chrono>
#include "RmgTimer.h"


std::chrono::high_resolution_clock::time_point start_time;
std::chrono::high_resolution_clock::time_point end_time;
const char *sname;
bool RmgTimer::enabled = true;
std::unordered_map<std::string, double> RmgTimer::timings[MAX_RMG_THREADS+1];



RmgTimer::RmgTimer(const char *what) 
{
    if(!RmgTimer::enabled) return;
    start_time = std::chrono::high_resolution_clock::now();
    sname = what;
}

RmgTimer::~RmgTimer(void) 
{
    if(!RmgTimer::enabled) return;
    if(!sname) return;  // Bit of a hack for the temporary C interfaces
    end_time = std::chrono::high_resolution_clock::now();

    int tid;
    BaseThread *T = BaseThread::getBaseThread(0);
    BaseThreadControl *Tptr;
    Tptr = T->get_thread_control();

    // if Tptr is null that means we are being called from the main program
    if(!Tptr) {
        tid = 0;
    }
    else {
        tid = T->get_thread_tid();
        if(tid < 0) {
            tid = 0;
        }
        else {
            tid++;      // Increment by 1 since we store the main thread values in the zeroth slot
        }
    }

    std::chrono::duration<double, std::milli> time_length = end_time - start_time;
    if(timings[tid].count(sname)) {
        RmgTimer::timings[tid][sname] += time_length.count() / 1000.0;
    }
    else {
        RmgTimer::timings[tid][sname] += time_length.count() / 1000.0;
    }
}

void RmgTimer::enable(void)
{
    RmgTimer::enabled = true;
}
void RmgTimer::disable(void)
{
    RmgTimer::enabled = false;
}

std::unordered_map<std::string, double> *RmgTimer::get_map(void)
{
    return RmgTimer::timings;
}

// Temporary until C to C++ migration is completed. Use carefully!
extern "C" void *BeginRmgTimer(const char *what) 
{
    if(!RmgTimer::enabled) return NULL;
    RmgTimer *RT = new RmgTimer(what);
    return (void *)RT;
}
extern "C" void EndRmgTimer(void *ptr) 
{
    if(!RmgTimer::enabled) return;
    RmgTimer *RT = (RmgTimer *)ptr;
    delete(RT);
}

double MyCrtc (void)
{
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> time_length = t1 - t1.min();
    return time_length.count() / 1000.0;

}

extern "C" double my_crtc (void)
{
    return MyCrtc();
}


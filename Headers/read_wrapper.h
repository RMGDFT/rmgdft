#pragma once
#ifndef read_wrapper_H
#define read_wrapper_H 1

#include <stdio.h>
#include <unistd.h>

static void read_double (int fhand, double * rp, int count);
static void read_int (int fhand, int *ip, int count);

static void read_double (int fhand, double * rp, int count)
{
    ssize_t wanted = sizeof (double) * (ssize_t)count;
    ssize_t size = read (fhand, rp, wanted);
    if(size != wanted)
    {
        std::cout << " wanted and readsize " << wanted <<" != " << size << std::endl;
        rmg::error("error reading");
    }


}

static void read_int (int fhand, int *ip, int count)
{
    int size = count * sizeof (int);
    if (size != read (fhand, ip, size))
        rmg::error("error reading");
}

#endif

#pragma once
#ifndef write_wrapper_H
#define write_wrapper_H 1

#include <stdio.h>
#include <unistd.h>

static size_t totalsize;
static void write_double (int fh, double * rp, int count);
static void write_int (int fh, int *ip, int count);

static void write_double (int fh, double * rp, int count)
{
    int size;

    size = count * sizeof (double);
    rmg::writefile(fh, rp, size);
    totalsize += size;
}


static void write_int (int fh, int *ip, int count)
{
    int size;

    size = count * sizeof (int);
    rmg::writefile(fh, ip, size);

    totalsize += size;
}
#endif

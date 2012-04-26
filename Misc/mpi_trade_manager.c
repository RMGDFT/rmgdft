#if 0
/*

  This code implements an mpi request manager for tradeimages. 




*/
#define _GNU_SOURCE

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <signal.h>
#include <sys/types.h>
#include <unistd.h>
#include <semaphore.h>
#include <errno.h>

#include "main.h"
#include "hybrid.h"

// Pthread ID of manager thread
static pthread_t manager_thread;
static pthread_cond_t mgr_cv;
static pthread_mutex_t mgr_mutex = PTHREAD_MUTEX_INITIALIZER;



// Rank of target node based on offsets from current node
int target_node[3][3][3];

#define RMG_MPI_RING_BUFFER_SIZE  (4 * THREADS_PER_NODE * 27)   // Should be large enough for max active MPI operations
#define RMG_MPI_SLOT_OPEN 1
#define RMG_MPI_SLOT_QUEUED 2
#define RMG_MPI_SLOT_WAITING 3
#define RMG_MPI_SLOT_COMPLETED 4

#define RMG_MPI_ISEND 1
#define RMG_MPI_IRECV 2

typedef struct {

  // RMG_MPI_SLOT_OPEN means no request in this slot.
  // RMG_MPI_SLOT_PENDING means request submitted to queue but transfer not initiated
  // RMG_MPI_SLOT_WAITING means request initiated and waiting completion
  int status;  // slot available or pending completion RMG_MPI_SLOT_OPEN,RMG_MPI_SLOT_PENDING

  MPI_Request request;
  MPI_Status mpi_status;
  int mpi_request_flag;

  REAL *buf;   // send or receive buffer
  int count;   // number of data points of type real
  int target;  // node rank of target
  int tid;     // Id of thread the request is associated with
  int tag;     // assigned by the manager thread
  int type;    // RMG_MPI_ISEND or RMG_MPI_IRECV

  MPI_Comm comm;
} RMG_MPI_trade_t;


// Each thread can have at most 26 send or receive ops outstanding (1,1,1) is an error
// since we map -1 to the first element in the array, 0 to the second and 1 to third. So
// a trade with (-1,-1,-1) is the [0][0][0] array element and (1,1,1) is [1][1][1] which is
// ourself.
RMG_MPI_trade_t RMG_MPI_RecvList[THREADS_PER_NODE][3][3][3];
RMG_MPI_trade_t RMG_MPI_SendList[THREADS_PER_NODE][3][3][3];


// Each thread uses these to wait
typedef struct {
  pthread_mutex_t lock;
  pthread_cond_t cv;
  volatile int recv_count;
  volatile int send_count;
} RMG_MPI_lock_t;

RMG_MPI_lock_t RMG_MPI_locks[THREADS_PER_NODE];


// Management structures for async MPI requests that have
// been queued but not yet executed.
static volatile int rb_head=0;       // Can be modified by the manager thread or the wfunc threads so both protect with mutex
static volatile int rb_tail=0;       // Only be modified by manager thread
static volatile int rb_exec_tail=0;  // Only modified by manager thread
static pthread_mutex_t rb_mutex = PTHREAD_MUTEX_INITIALIZER;
RMG_MPI_trade_t *RMG_MPI_RingBuffer[8 * THREADS_PER_NODE * 27];


void RMG_MPI_push_req(RMG_MPI_trade_t *req) {
   pthread_mutex_lock(&rb_mutex);
   //printf("Manager: pushing %d %d %d %d\n",RMG_MPI_RING_BUFFER_SIZE, rb_head,rb_tail,rb_exec_tail);
   //fflush(NULL);
   RMG_MPI_RingBuffer[rb_head] = req;
   rb_head++;
   rb_head = rb_head % RMG_MPI_RING_BUFFER_SIZE;
   pthread_mutex_unlock(&rb_mutex);
}


RMG_MPI_trade_t *RMG_MPI_pull_req(void) {
   RMG_MPI_trade_t *reqptr;

   
   pthread_mutex_lock(&rb_mutex);
   //printf("Manager: pulling %d %d %d\n",rb_head,rb_tail,rb_exec_tail);
   //fflush(NULL);
   if(rb_head == rb_tail) {
       pthread_mutex_unlock(&rb_mutex);
       return NULL;
   }
   reqptr = RMG_MPI_RingBuffer[rb_tail];
   rb_tail++;
   rb_tail = rb_tail % RMG_MPI_RING_BUFFER_SIZE;
   pthread_mutex_unlock(&rb_mutex);

   pthread_mutex_lock(&RMG_MPI_locks[reqptr->tid].lock);
   if(reqptr->type == RMG_MPI_IRECV) {
       MPI_Irecv(reqptr->buf, reqptr->count, MPI_DOUBLE, reqptr->target,
                  reqptr->tag, reqptr->comm, &reqptr->request);
       reqptr->status = RMG_MPI_SLOT_WAITING;
   }
   else {
       MPI_Isend(reqptr->buf, reqptr->count, MPI_DOUBLE, reqptr->target,
                      reqptr->tag, reqptr->comm, &reqptr->request);
       reqptr->status = RMG_MPI_SLOT_WAITING;
   }
   pthread_mutex_unlock(&RMG_MPI_locks[reqptr->tid].lock);

   // Recursive call
   RMG_MPI_pull_req();
   return reqptr;
}


// This function is used to insert the request into the queues.
void RMG_MPI_trade(REAL *buf, int count, int type, int pe_x_offset, int pe_y_offset, int pe_z_offset, MPI_Comm comm, int tag)
{
    int tid;
    RMG_MPI_trade_t *reqptr; 


    // Incement offsets so they can act as array indices into send and recv lists
    pe_x_offset++;
    pe_y_offset++;
    pe_z_offset++;

    // Tag is based on tid in the lower 8 bits with the array

    tid = get_thread_tid();
    if(tid == -1) tid = 0;

    pthread_mutex_lock(&RMG_MPI_locks[tid].lock);

    if(type == RMG_MPI_ISEND) {
        reqptr = &RMG_MPI_SendList[tid][pe_x_offset][pe_y_offset][pe_z_offset];
        RMG_MPI_locks[tid].send_count++;
        // Check if the slot is open. If not then it is an error so abort
        if(reqptr->status != RMG_MPI_SLOT_OPEN) {
            printf("Error in RMG_MPI_trade (send). Slot already occupied. tid=%d\n", tid);
            fflush(NULL);
            exit(0);
        }

    }
    else {
        reqptr = &RMG_MPI_RecvList[tid][pe_x_offset][pe_y_offset][pe_z_offset];
        RMG_MPI_locks[tid].recv_count++;
        // Check if the slot is open. If not then it is an error so abort
        if(reqptr->status != RMG_MPI_SLOT_OPEN) {
            printf("Error in RMG_MPI_Isend (recv). Slot already occupied. tid=%d\n", tid);
            fflush(NULL);
            exit(0);
        }

    }

    reqptr->tid = tid;
    reqptr->tag = (tag<<16) + tid;
    reqptr->buf = buf;
    reqptr->count = count;
    reqptr->target = target_node[pe_x_offset][pe_y_offset][pe_z_offset];
    reqptr->status = RMG_MPI_SLOT_QUEUED;
    reqptr->comm = comm;
    reqptr->type = type;
    RMG_MPI_push_req(reqptr);
    pthread_mutex_unlock(&RMG_MPI_locks[tid].lock);

    pthread_cond_signal(&mgr_cv);
    sched_yield();

}


// Waits for all pending requests for this thread to complete
void RMG_MPI_trade_waitall(void) 
{
    int tid, rc;
    struct timespec   ts;
    struct timeval   tp;

    rc =  gettimeofday(&tp, NULL);
    ts.tv_sec  = tp.tv_sec;
    ts.tv_nsec = tp.tv_usec * 1000 + 50000;
//    ts.tv_sec += 2;
 
    tid = get_thread_tid();
    if(tid == -1) tid = 0;

    // First check if any requests are pending and if no then just return
    pthread_mutex_lock(&RMG_MPI_locks[tid].lock);
    if((RMG_MPI_locks[tid].send_count == 0) && (RMG_MPI_locks[tid].recv_count == 0)) {
        pthread_mutex_unlock(&RMG_MPI_locks[tid].lock);
        return;
    }

    // Else wait wake up the main thread and wait on the condition variable
//    pthread_cond_timedwait(&RMG_MPI_locks[tid].cv, &RMG_MPI_locks[tid].lock, &ts);
    pthread_cond_signal(&mgr_cv);
    pthread_cond_wait(&RMG_MPI_locks[tid].cv, &RMG_MPI_locks[tid].lock);
    pthread_mutex_unlock(&RMG_MPI_locks[tid].lock);
}


// Waits for all pending recv requests for this thread to complete
void RMG_MPI_trade_recvwait(void) 
{
    int tid;

    tid = get_thread_tid();
    if(tid == -1) tid = 0;

    // First check if any requests are pending and if no then just return
    pthread_mutex_lock(&RMG_MPI_locks[tid].lock);
    if(RMG_MPI_locks[tid].recv_count == 0) {
        pthread_mutex_unlock(&RMG_MPI_locks[tid].lock);
        return;
    }

    // Else wait on the condition variable
    pthread_cond_wait(&RMG_MPI_locks[tid].cv, &RMG_MPI_locks[tid].lock);
    pthread_mutex_unlock(&RMG_MPI_locks[tid].lock);
}



//
void trade_images_manager(void *s)
{
    int thread_idx, req_idx, rc, retval;
    RMG_MPI_trade_t *reqptr, *reqptr1;
    struct timespec rqtp;
    int qidx, qidx1, qlen;
    struct timespec   ts;
    struct timeval   tp;

    int ix, iy, iz;
    int pe_x, pe_y, pe_z;
    int t_pe_x, t_pe_y, t_pe_z;

    rqtp.tv_sec = 0;
    rqtp.tv_nsec = 5000;

    // Set up the target node array
    pe2xyz(pct.gridpe, &pe_x, &pe_y, &pe_z);
    for(ix = -1;ix <= 1;ix++) {
        t_pe_x = pe_x + ix;
        if(t_pe_x < 0) t_pe_x = PE_X - 1;
        if(t_pe_x == PE_X) t_pe_x = 0;
        for(iy = -1;iy <= 1;iy++) {
            t_pe_y = pe_y + iy;
            if(t_pe_y < 0) t_pe_y = PE_Y - 1;
            if(t_pe_y == PE_Y) t_pe_y = 0;

            for(iz = -1;iz <= 1;iz++) {
                t_pe_z = pe_z + iz;
                if(t_pe_z < 0) t_pe_z = PE_Z - 1;
                if(t_pe_z == PE_Z) t_pe_z = 0;
                target_node[ix+1][iy+1][iz+1] = t_pe_x*PE_Y*PE_Z + t_pe_y*PE_Z + t_pe_z;
            }
        }
    } // end for


    // Initialize the lock and trade structures
    for(thread_idx = 0;thread_idx < THREADS_PER_NODE;thread_idx++) {
        RMG_MPI_locks[thread_idx].recv_count = 0;
        RMG_MPI_locks[thread_idx].send_count = 0;
        pthread_mutex_init(&RMG_MPI_locks[thread_idx].lock, NULL);
        pthread_cond_init(&RMG_MPI_locks[thread_idx].cv, NULL);
        reqptr = &RMG_MPI_SendList[thread_idx][0][0][0];
        for(req_idx = 0;req_idx < 27;req_idx++) {
            reqptr->status = RMG_MPI_SLOT_OPEN;
            reqptr++;
        }
        reqptr = &RMG_MPI_RecvList[thread_idx][0][0][0];
        for(req_idx = 0;req_idx < 27;req_idx++) {
            reqptr->status = RMG_MPI_SLOT_OPEN;
            reqptr++;
        }
    } // end for

    // Initialize our condition variable
    pthread_cond_init(&mgr_cv, NULL);

    // Set our id in the global so threads have access to it
    manager_thread = pthread_self();

    // Loop forever
    while(1) {

        RMG_MPI_pull_req();

//        nanosleep(&rqtp, NULL);
    rc =  gettimeofday(&tp, NULL);
    ts.tv_sec  = tp.tv_sec;
    ts.tv_nsec = tp.tv_usec * 1000 + 5000;
//    ts.tv_sec += 1;
        pthread_mutex_lock(&mgr_mutex);
        retval = pthread_cond_timedwait(&mgr_cv, &mgr_mutex, &ts);
        pthread_mutex_unlock(&mgr_mutex);


        pthread_mutex_lock(&rb_mutex);
        qlen = rb_tail - rb_exec_tail;
        if(qlen < 0) qlen += RMG_MPI_RING_BUFFER_SIZE;
        qidx1 = rb_exec_tail;
        pthread_mutex_unlock(&rb_mutex);
        for(qidx = 0;qidx < qlen;qidx++) {

        RMG_MPI_pull_req();
            //printf("Manager: looping %d %d %d %d %d\n",rb_head,rb_tail,rb_exec_tail,qlen,qidx1);
            //fflush(NULL);
            reqptr = RMG_MPI_RingBuffer[qidx1];
            pthread_mutex_lock(&RMG_MPI_locks[reqptr->tid].lock);
            if(reqptr->status == RMG_MPI_SLOT_WAITING) {
                MPI_Test(&reqptr->request, &reqptr->mpi_request_flag, &reqptr->mpi_status);
                if(reqptr->mpi_request_flag) {
                    //printf("Manager: tested head=%d tail=%d exec_tail=%d qlen=%d qidx=%d qidx1=%d\n",rb_head,rb_tail,rb_exec_tail,qlen,qidx,qidx1);
                    //fflush(NULL);
                    reqptr->status = RMG_MPI_SLOT_OPEN;
                    if(reqptr->type == RMG_MPI_IRECV) {
                        RMG_MPI_locks[reqptr->tid].recv_count--;
                    }
                    else {
                        RMG_MPI_locks[reqptr->tid].send_count--;
                    }
                    if(qidx != 0) {
                        
                        //printf("Manager: swapping head=%d tail=%d exec_tail=%d qlen=%d qidx=%d qidx1=%d\n",rb_head,rb_tail,rb_exec_tail,qlen,qidx,qidx1);
                        //fflush(NULL);
                        RMG_MPI_RingBuffer[qidx1] = RMG_MPI_RingBuffer[rb_exec_tail];

                    }
                    rb_exec_tail++;
                    rb_exec_tail = rb_exec_tail % RMG_MPI_RING_BUFFER_SIZE;
                    if((RMG_MPI_locks[reqptr->tid].recv_count == 0) && (RMG_MPI_locks[reqptr->tid].send_count == 0)) {
                        pthread_cond_signal(&RMG_MPI_locks[reqptr->tid].cv);
                    }
                }
            }
            pthread_mutex_unlock(&RMG_MPI_locks[reqptr->tid].lock);
            qidx1++;
            qidx1 = qidx1 % RMG_MPI_RING_BUFFER_SIZE;
        }

    }
}
#endif

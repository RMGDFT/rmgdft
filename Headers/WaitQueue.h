#ifndef RMG_WaitQueue_H
#define RMG_WaitQueue_H 1

#include <queue>
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>

template<typename T>
class WaitQueue
{
private:
    std::queue<T> q;
    mutable boost::mutex m;
    boost::condition_variable cv;
public:
    void push(T const& data)
    {
        boost::mutex::scoped_lock lock(m);
        q.push(data);
        lock.unlock();
        cv.notify_one();
    }

    bool is_empty() const
    {
        boost::mutex::scoped_lock lock(m);
        return q.empty();
    }

    bool try_pop(T& popped_value)
    {
        boost::mutex::scoped_lock lock(m);
        if(q.is_empty())
        {
            return false;
        }
        
        popped_value=q.front();
        q.pop();
        return true;
    }

    void wait_and_pop(T& popped_value)
    {
        boost::mutex::scoped_lock lock(m);
        while(q.empty())
        {
            cv.wait(lock);
        }
        
        popped_value=q.front();
        q.pop();
    }

};

#endif

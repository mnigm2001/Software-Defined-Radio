#include <iostream>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>

//For random number generation
#define RANDOM_BOUND 1000
#define RANDOM_SEED 0

// a define and a variable used to keep track when to terminate the threads
#define TOTAL_ELEMS 10
int generated_elems = 0; // global variable (only for illustrative purposes)

// the size of the queue
#define QUEUE_ELEMS 3

/*
struct TSYNC {
  std::mutex *mutex;
  std::condition_variable *cvar;
};*/


// some generaton function (for illustration )
void my_generate(int& elem)
{
  elem = rand() % (RANDOM_BOUND);
  elem -= int(RANDOM_BOUND/2);
  std::cerr << "Generated " << elem << std::endl;
}

// some compute function
void my_compute(int& elem)
{
  std::cerr << "Recieved " << elem << std::endl;
}

void my_produce(std::queue<int> &my_queue, std::mutex& my_mutex, \
  std::condition_variable& my_cvar)
{
  // thread runs until a custom end cond occurs
  while(generated_elems < TOTAL_ELEMS)
  {
    //generation and custom processing of elements before
    //passing them to another thread through a queue
    int elem;
    //any processing idependent of other threads should not use locks
    my_generate(elem);
    //bookkeeping for the end cond (very custom)
    generated_elems++;

    //pass data to the queue using mutexes and
    // condition variables as a sync mechanism
    std::unique_lock<std::mutex> my_lock(my_mutex);
    //check if the queue is FULL - needed to be done in a custom way
    while(my_queue.size() >= QUEUE_ELEMS){
      my_cvar.wait(my_lock);
    }

    my_queue.push(elem);
    my_cvar.notify_one();
    my_lock.unlock();
  }
}

//void my_consume(std::queue<int> &my_queue, std::mutex& my_mutex, std::condition_variable& my_cvar)
void my_consume(std::queue<int> &my_queue)
{
  
  //threads ofter operate in infinite loops and
  // they break when a custom end cond occurs
  while(1)
  {
    // extract data from the queue using mutexes and
    // cond variables as a sync mechanism
    std::unique_lock<std::mutex> my_lock(my_mutex);
    while(my_queue.empty()){
      my_cvar.wait(my_lock);
    }

    int elem = my_queue.front();
    my_queue.pop();
    my_cvar.notify_one();
    my_lock.unlock();

    // use the extracted data from the queue to do some custom processing
    // the mutex should be unclocked to enable parallel execution on multicore
    my_compute(elem);
    // check the end cond (custom) to terminate the thread
    if ((generated_elems == TOTAL_ELEMS) && my_queue.empty()){
      break;
    }
  }
}


int main() {
  // initialize the random generator (illustrative)
  srand(RANDOM_SEED);

  std::cout << "Hardware support for " << std::thread:hardware_concurrency() << "parallel threads" << std::endl;


  // context used for thread synch through queue
  std::queue<int> my_queue;
  std::mutex my_mutex;
  std::condition_variable my_cvar;

  /*
  struct TSYNC T1_SYNC = {
    &my_mutex,
    &my_cvar
  };*/

  // producer thread
  std::thread tp = std::thread(my_produce, std::ref(my_queue), \
    std::ref(my_mutex), std::ref(my_cvar));

  /*
  // producer thread
  std::thread tp = std::thread(my_produce, std::ref(my_queue), \
    std::ref(T1_SYNC));*/

  // consumer thread
  std::thread tc = std::thread(my_consume, std::ref(my_queue));


  /*
  // consumer thread
  std::thread tc = std::thread(my_consume, std::ref(my_queue), \
    std::ref(my_mutex), std::ref(my_cvar));*/

  tp.join();
  tc.join();

  return 0;
}

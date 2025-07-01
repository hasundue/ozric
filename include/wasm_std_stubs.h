// WASM Standard Library Stubs
// This header provides minimal implementations of standard library features
// that aren't available in WASM but are needed by Ceres headers
// Include this BEFORE any standard library headers

#pragma once

#ifdef __wasm__

// Prevent standard library headers from being included
#define _MUTEX_
#define _THREAD_
#define _CONDITION_VARIABLE_
#define _ATOMIC_

// Provide our stub implementations
namespace std {
    // Mutex stub - no-op for single-threaded WASM
    class mutex {
    public:
        void lock() {}
        void unlock() {}
        bool try_lock() { return true; }
    };
    
    // Lock guard stub 
    template<typename Mutex>
    class lock_guard {
    public:
        explicit lock_guard(Mutex&) {}
        ~lock_guard() {}
    };
    
    // Condition variable stub
    class condition_variable {
    public:
        void notify_one() {}
        void notify_all() {}
        template<typename Lock>
        void wait(Lock&) {}
        template<typename Lock, typename Predicate>
        void wait(Lock&, Predicate) {}
    };
    
    // Thread stub 
    class thread {
    public:
        template<typename Function, typename... Args>
        thread(Function&&, Args&&...) {}
        
        void join() {}
        void detach() {}
        bool joinable() const { return false; }
        
        class id {};
        id get_id() const { return id{}; }
        
        static unsigned hardware_concurrency() { return 1; }
    };
    
    // Atomic stub - single-threaded, no synchronization needed
    template<typename T>
    class atomic {
        T value_;
    public:
        atomic() = default;
        atomic(T desired) : value_(desired) {}
        
        T load() const { return value_; }
        void store(T desired) { value_ = desired; }
        T exchange(T desired) { 
            T old = value_; 
            value_ = desired; 
            return old; 
        }
        
        operator T() const { return load(); }
        T operator=(T desired) { store(desired); return desired; }
    };

    // Memory order enum stub
    enum class memory_order {
        relaxed, consume, acquire, release, acq_rel, seq_cst
    };
}

#endif // __wasm__
# Ceres-solver Code Analysis

This document contains findings from analyzing the Ceres-solver source code integrated into this project.

## Source Code Location

The Ceres-solver 2.2.0 source code is located in Zig's dependency cache:
- **Path**: `/home/hasundue/.cache/zig/p/N-V-__8AAAwZyQAma1OByLx2AKj3upoY7FxoS7bQ1LBLR2m3/`
- **Hash**: `N-V-__8AAAwZyQAma1OByLx2AKj3upoY7FxoS7bQ1LBLR2m3` (from build.zig.zon)

## Parallel CostFunction::Evaluate Calls

**Evidence of threading**: Ceres calls `CostFunction::Evaluate` in parallel with threads.

**Key files**:
- `/internal/ceres/program_evaluator.h:186-300` - Core parallel evaluation loop
- `/internal/ceres/parallel_for.h` - Distributes work across threads  
- `/internal/ceres/thread_pool.cc` - Work-stealing thread pool
- `/internal/ceres/residual_block.cc:112` - Where `CostFunction::Evaluate` is called

**Call chain**: `ProgramEvaluator::Evaluate()` → `ParallelFor()` → `ResidualBlock::Evaluate()` → `CostFunction::Evaluate()`

**Threading architecture**:
- Work-stealing thread pool with shared task queues
- Each residual block evaluated independently across threads
- Thread-local scratch space prevents contention
- Uses `std::thread::hardware_concurrency()` for optimal thread count
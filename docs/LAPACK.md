# LAPACK Integration and Implementation

This document covers various aspects of LAPACK (Linear Algebra Package) integration and implementation in the Ozric project.

## Overview

LAPACK is a library of Fortran 90 subroutines for solving systems of linear equations, linear least squares problems, eigenvalue problems, and singular value decomposition. Our project implements selected LAPACK-compatible routines in Zig for numerical optimization in density functional theory.

## Project Integration

### Current Implementation Status

- **Implemented**: `dsbmv` - Symmetric banded matrix-vector multiplication
- **Planned**: Additional BLAS Level 2 and Level 3 routines as needed
- **Target**: LAPACK-compatible API with Zig's memory safety and performance

### Design Principles

1. **API Compatibility**: Function signatures match LAPACK conventions
2. **Memory Safety**: Zig's compile-time and runtime safety features
3. **Performance**: Zero-cost abstractions and SIMD optimization potential
4. **Modularity**: Individual routines can be used independently

## LAPACK Naming Conventions

LAPACK uses systematic naming conventions for its routines:

### Routine Name Structure: `[precision][matrix_type][operation]`

- **Precision**: `S` (single), `D` (double), `C` (complex), `Z` (double complex)
- **Matrix Type**: `GE` (general), `SY` (symmetric), `HE` (hermitian), `GB` (banded), `SB` (symmetric banded)
- **Operation**: `MV` (matrix-vector), `MM` (matrix-matrix), `SV` (solve), etc.

### Examples
- `dsbmv`: **D**ouble precision **S**ymmetric **B**anded **M**atrix-**V**ector multiplication
- `dgemv`: **D**ouble precision **GE**neral **M**atrix-**V**ector multiplication  
- `dsyev`: **D**ouble precision **SY**mmetric **E**igenvalue decomposition

## Matrix Storage Formats

LAPACK uses various storage formats optimized for different matrix types:

### 1. Full Storage (General Matrices)
- Standard row-major or column-major storage
- Used for general dense matrices
- Simple indexing: `A[i + j * lda]`

### 2. Packed Storage (Symmetric/Triangular)
- Stores only upper or lower triangle
- Saves ~50% memory for symmetric matrices
- More complex indexing required

### 3. Banded Storage
- For matrices with non-zero elements only near the diagonal
- Stores only the band containing non-zero elements
- Significant memory savings for sparse banded matrices

## Symmetric Banded Matrix Storage (dsbmv)

### Storage Format

For a symmetric banded matrix with bandwidth `k`, the storage uses:
- **Rows**: Represent diagonals (from superdiagonal k down to main diagonal)
- **Columns**: Represent matrix columns
- **Leading dimension**: `lda = k + 1`

### Indexing Formula

For upper triangle storage, matrix element `A(i,j)` is stored at:
```
a[k + i - j + j * lda]
```

### Example: 3×3 Matrix with Bandwidth 1

Consider the symmetric tridiagonal matrix:
```
A = [[2, 1, 0],
     [1, 2, 1],
     [0, 1, 2]]
```

With `k = 1` and `lda = 2`, the storage layout becomes:

```
Storage Array Layout (column-major):
                     Col 0  Col 1  Col 2
Row 0 (superdiag):  [unused, A(0,1), A(1,2)]  = [0, 1, 1]
Row 1 (diagonal):   [A(0,0), A(1,1), A(2,2)] = [2, 2, 2]
```

Resulting in the array: `[0, 2, 1, 2, 1, 2]`
(Reading column-major: unused, A(0,0), A(0,1), A(1,1), A(1,2), A(2,2))

### Element Mapping

Using the indexing formula `a[k + i - j + j * lda]` with `k=1`, `lda=2`:
- `A(0,0)` → `a[1 + 0 - 0 + 0*2]` = `a[1]` = 2
- `A(0,1)` → `a[1 + 0 - 1 + 1*2]` = `a[2]` = 1  
- `A(1,1)` → `a[1 + 1 - 1 + 1*2]` = `a[3]` = 2
- `A(1,2)` → `a[1 + 1 - 2 + 2*2]` = `a[4]` = 1
- `A(2,2)` → `a[1 + 2 - 2 + 2*2]` = `a[5]` = 2

### Why is the First Element Unused?

The first element `a[0]` would correspond to `A(0,-1)` using the indexing formula, which doesn't exist in the matrix. This is a byproduct of having a consistent, simple indexing formula that works for all valid elements in the banded region.

## Parameter Conventions

### Common Parameters

- **`uplo`**: 'U' (upper triangle) or 'L' (lower triangle)
- **`n`**: Matrix dimension
- **`k`**: Number of super/sub-diagonals in banded matrices
- **`lda`**: Leading dimension of array
- **`incx`, `incy`**: Stride between elements in vectors
- **`alpha`, `beta`**: Scalar coefficients

### Our Zig Implementation

```zig
pub const UPLO = enum {
    U, // Upper triangle storage
    L, // Lower triangle storage (not implemented)
};

pub fn dsbmv(
    comptime uplo: UPLO,
    n: usize,
    k: usize,
    alpha: f64,
    a: []const f64,
    lda: usize,
    x: []const f64,
    incx: isize,
    beta: f64,
    y: []f64,
    incy: isize,
) void
```

## Error Handling

### LAPACK Approach
- Uses `info` parameter to return error codes
- Positive values indicate specific errors
- Zero indicates successful completion

### Our Zig Approach
- **Compile-time errors**: For unsupported operations (e.g., lower triangle storage)
- **Runtime safety**: Bounds checking and memory safety
- **Explicit error types**: Using Zig's error system where appropriate

## Performance Considerations

### Memory Access Patterns
- **Column-major**: LAPACK's default for Fortran compatibility
- **Cache efficiency**: Optimized for typical linear algebra operations
- **SIMD potential**: Our implementation includes SIMD optimization framework

### Optimization Opportunities
1. **Vectorization**: SIMD instructions for parallel operations
2. **Blocking**: Cache-friendly algorithms for large matrices
3. **Specialization**: Compile-time optimization for specific matrix sizes

## Testing and Validation

### Test Coverage
- **Basic functionality**: Standard matrix operations
- **Edge cases**: Single elements, zero matrices, special parameters
- **Parameter combinations**: Different alpha/beta values, strides
- **Numerical precision**: Small and large number handling

### LAPACK Compatibility
- **Reference implementations**: Compare against standard LAPACK
- **Numerical accuracy**: Maintain precision standards
- **API consistency**: Match expected behavior exactly

## Future Extensions

### Planned Implementations
- **`dgemv`**: General matrix-vector multiplication
- **`dsymv`**: Symmetric matrix-vector multiplication
- **`dger`**: General rank-1 update
- **`dsyr`**: Symmetric rank-1 update

### Integration with DFT
- **Density matrices**: Symmetric matrix operations
- **Hamiltonian operations**: Banded matrix structures common in DFT
- **Eigenvalue problems**: Future eigensolvers for electronic structure

## References

- [LAPACK Users' Guide](http://www.netlib.org/lapack/lug/)
- [LAPACK Reference Manual](http://www.netlib.org/lapack/explore-html/)
- [BLAS Technical Documentation](http://www.netlib.org/blas/)
- [Matrix Storage Formats in LAPACK](http://www.netlib.org/lapack/lug/node121.html)
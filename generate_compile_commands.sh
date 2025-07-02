#!/usr/bin/env bash

# Generate compile_commands.json for clangd
# This script creates a compile database for C++ files in the project

set -e

echo "Generating compile_commands.json for clangd..."

# Get absolute project path
PROJECT_DIR="$(pwd)"

# Find ceres dependency path
CERES_PATH=$(find ~/.cache/zig/p -name "N-V-__8AAAwZyQAma1OByLx2AKj3upoY7FxoS7bQ1LBLR2m3" -type d 2>/dev/null | head -1)
EIGEN_PATH=$(find ~/.cache/zig/p -name "N-V-__8AABPI6wAbBLaww8cKO-HchZ0ft-BxUOEOyQtTSccD" -type d 2>/dev/null | head -1)

if [ -z "$CERES_PATH" ] || [ -z "$EIGEN_PATH" ]; then
    echo "Error: Ceres or Eigen dependencies not found. Run 'zig build' first to fetch dependencies."
    exit 1
fi

# Create compile_commands.json with proper JSON escaping
cat > compile_commands.json << 'EOF'
[
  {
    "directory": "PROJECT_DIR_PLACEHOLDER",
    "command": "clang++ -std=c++17 -DCERES_NO_SUITESPARSE -DCERES_NO_CXSPARSE -DCERES_EXPORT= -DCERES_NO_EXPORT= -DCERES_NO_PROTOCOL_BUFFERS -DCERES_NO_THREADS -DCERES_NO_CUDA -DCERES_NO_LAPACK -DCERES_METIS_VERSION=\"5.1.0\" -DMINIGLOG -DCERES_RESTRICT_SCHUR_SPECIALIZATION -DCERES_NO_ACCELERATE_SPARSE -IPROJECT_DIR_PLACEHOLDER/include -ICERES_PATH_PLACEHOLDER/include -ICERES_PATH_PLACEHOLDER -ICERES_PATH_PLACEHOLDER/internal -ICERES_PATH_PLACEHOLDER/config -ICERES_PATH_PLACEHOLDER/internal/ceres/miniglog -IEIGEN_PATH_PLACEHOLDER -c src/solver.cc",
    "file": "src/solver.cc"
  }
]
EOF

# Replace placeholders with actual paths
sed -i "s|PROJECT_DIR_PLACEHOLDER|${PROJECT_DIR}|g" compile_commands.json
sed -i "s|CERES_PATH_PLACEHOLDER|${CERES_PATH}|g" compile_commands.json  
sed -i "s|EIGEN_PATH_PLACEHOLDER|${EIGEN_PATH}|g" compile_commands.json

echo "compile_commands.json generated successfully!"
echo "You can now use clangd for C++ language support."
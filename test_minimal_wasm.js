#!/usr/bin/env node

console.log('🚀 Testing minimal WASM build for comparison...');

const { execSync } = require('child_process');

try {
    // Test the minimal WASM with wasmtime
    console.log('⏳ Running minimal WASM with wasmtime...');
    const result = execSync('nix develop --command wasmtime zig-out/bin/ozric_wasm.wasm', 
                           { encoding: 'utf8', timeout: 5000 });
    console.log('✅ Minimal WASM result:', result);
} catch (error) {
    console.log('❌ Minimal WASM error:', error.message);
}

console.log('\n📊 File size comparison:');
console.log('Minimal WASM:', execSync('ls -lh zig-out/bin/ozric_wasm.wasm').toString().split(/\s+/)[4]);
console.log('Threading WASM:', execSync('ls -lh zig-out/bin/ozric_wasm_threads.wasm').toString().split(/\s+/)[4]);
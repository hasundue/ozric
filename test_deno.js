#!/usr/bin/env deno run --allow-read

console.log('🚀 Testing Ceres WASM with Deno...');

try {
    // Load the WASM file directly
    const wasmBytes = await Deno.readFile('./zig-out/bin/ozric_wasm_threads.wasm');
    console.log(`📦 WASM file loaded: ${wasmBytes.length} bytes`);
    
    // Try to instantiate it
    const wasmModule = await WebAssembly.instantiate(wasmBytes, {
        env: {},
        wasi_snapshot_preview1: {}
    });
    
    console.log('✅ WASM module instantiated');
    console.log('📋 Exported functions:', Object.keys(wasmModule.instance.exports));
    
} catch (error) {
    console.error('❌ Error:', error.message);
}
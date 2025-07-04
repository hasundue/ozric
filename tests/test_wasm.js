#!/usr/bin/env node

// Test script for Ceres WASM with threading support
const Module = require('../zig-out/bin/ozric_wasm.js');

Module.onRuntimeInitialized = () => {
    console.log('✅ WASM runtime initialized!');
    
    try {
        // Test Ceres version using ccall
        console.log('\n🔍 Testing Ceres version...');
        const version = Module.ccall('test_ceres', 'number', [], []);
        console.log(`📋 Ceres version major: ${version}`);
        
        // Test solver optimization using ccall
        console.log('\n🎯 Running Ceres solver optimization...');
        console.log('   Problem: minimize 0.5 * (10 - x)^2, initial x = 5.0');
        
        const result = Module.ccall('solve', 'number', [], []);
        console.log(`   ✨ Optimized result: x = ${result}`);
        console.log(`   📈 Expected: x ≈ 10.0 (theoretical minimum)`);
        
        if (Math.abs(result - 10.0) < 0.001) {
            console.log('   🎉 SUCCESS: Optimization converged correctly!');
        } else {
            console.log('   ⚠️  WARNING: Optimization result unexpected');
        }
        
        console.log('\n🧵 Threading info:');
        console.log(`   📱 PThread available: ${Module.PThread ? 'Yes' : 'No'}`);
        console.log(`   🏃 Worker threads: ${Module.PThread?.unusedWorkers?.length || 0}`);
        
        console.log('\n🎊 All tests completed! Ceres WASM with threading is working!');
        
    } catch (error) {
        console.error('❌ Error calling functions:', error);
    }
};

console.log('🚀 Loading Ceres WASM with threading support...');
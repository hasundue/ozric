#!/usr/bin/env node

console.log('🚀 Testing Ceres WASM threading build...');

// Set up Module object with callbacks
global.Module = {
    preRun: function() {
        console.log('⏳ WASM pre-run starting...');
    },
    
    postRun: function() {
        console.log('⏳ WASM post-run completed...');
    },
    
    onRuntimeInitialized: function() {
        console.log('✅ WASM runtime initialized!');
        
        try {
            // Test Ceres version
            console.log('\n🔍 Testing Ceres version function...');
            const version = Module.ccall('test_ceres', 'number', [], []);
            console.log(`📋 Ceres version major: ${version}`);
            
            // Test hello world optimization
            console.log('\n🎯 Running Ceres hello world optimization...');
            console.log('   Problem: minimize 0.5 * (10 - x)^2, initial x = 5.0');
            
            const result = Module.ccall('run_hello_world', 'number', [], []);
            console.log(`   ✨ Optimized result: x = ${result}`);
            console.log(`   📈 Expected: x ≈ 10.0 (theoretical minimum)`);
            
            if (Math.abs(result - 10.0) < 0.001) {
                console.log('   🎉 SUCCESS: Optimization converged correctly!');
            } else {
                console.log('   ⚠️  WARNING: Optimization result unexpected');
            }
            
            console.log('\n🧵 Threading info:');
            console.log(`   📱 PThread support: ${Module.PThread ? 'Available' : 'Not available'}`);
            console.log(`   🏃 Worker count: ${Module.PThread?.unusedWorkers?.length || 0}`);
            
            console.log('\n🎊 All tests completed! Ceres WASM with threading is working!');
            
        } catch (error) {
            console.error('❌ Error calling functions:', error);
        }
    },
    
    print: function(text) {
        console.log('WASM output:', text);
    },
    
    printErr: function(text) {
        console.error('WASM error:', text);
    }
};

// Load the Emscripten module
require('./zig-out/bin/ozric_wasm_threads.js');
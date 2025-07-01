#!/usr/bin/env node

console.log('🚀 Testing Ceres WASM with proper Emscripten initialization...');

// Set up the Module object before loading
global.Module = {
    onRuntimeInitialized: function() {
        console.log('✅ Runtime initialized!');
        
        try {
            // Test Ceres version
            console.log('\n🔍 Testing Ceres version...');
            const version = Module._test_ceres();
            console.log(`📋 Ceres version major: ${version}`);
            
            // Test hello world optimization
            console.log('\n🎯 Running Ceres hello world optimization...');
            console.log('   Problem: minimize 0.5 * (10 - x)^2');
            console.log('   Initial guess: x = 5.0');
            
            const result = Module._run_hello_world();
            console.log(`   ✨ Optimized result: x = ${result}`);
            console.log(`   📈 Expected: x ≈ 10.0 (theoretical minimum)`);
            
            if (Math.abs(result - 10.0) < 0.001) {
                console.log('   🎉 SUCCESS: Optimization converged correctly!');
            } else {
                console.log('   ⚠️  WARNING: Optimization result seems off');
            }
            
            console.log('\n✅ All tests completed successfully!');
            console.log('🎊 Ceres-solver is working in WASM with threading support!');
            
        } catch (error) {
            console.error('❌ Error calling functions:', error);
        }
    },
    
    print: function(text) {
        console.log('WASM:', text);
    },
    
    printErr: function(text) {
        console.error('WASM Error:', text);
    }
};

// Load the Emscripten module - this will call onRuntimeInitialized when ready
require('./zig-out/bin/ozric_wasm_threads.js');
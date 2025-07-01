#!/usr/bin/env node

console.log('🚀 Testing WASM threading build (synchronous approach)...');

// Disable threading for Node.js testing
global.Module = {
    noInitialRun: true,
    print: (text) => console.log('WASM:', text),
    printErr: (text) => console.error('WASM Error:', text),
    
    onRuntimeInitialized: function() {
        console.log('✅ Runtime initialized, testing functions...');
        
        try {
            // Check if functions exist
            if (typeof Module._test_ceres === 'function') {
                console.log('📋 test_ceres function found');
                const version = Module._test_ceres();
                console.log(`🔍 Ceres version: ${version}`);
            } else {
                console.log('❌ test_ceres function not found');
                console.log('Available functions:', Object.keys(Module).filter(k => k.startsWith('_')));
            }
            
            if (typeof Module._run_hello_world === 'function') {
                console.log('🎯 run_hello_world function found');
                console.log('⏳ Running optimization...');
                const result = Module._run_hello_world();
                console.log(`✨ Optimization result: x = ${result}`);
                
                if (Math.abs(result - 10.0) < 0.001) {
                    console.log('🎉 SUCCESS: Converged to correct solution!');
                } else {
                    console.log('⚠️ Unexpected result');
                }
            } else {
                console.log('❌ run_hello_world function not found');
            }
            
        } catch (error) {
            console.error('❌ Error calling functions:', error);
        }
    }
};

// Load the module
require('./zig-out/bin/ozric_wasm_threads.js');
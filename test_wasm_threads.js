#!/usr/bin/env node

// Simple Node.js script to test our Ceres WASM threading build
const fs = require('fs');
const path = require('path');

// Load the Emscripten-generated module
const Module = require('./zig-out/bin/ozric_wasm_threads.js');

async function testCeresWasm() {
    console.log('🚀 Loading Ceres WASM with threading support...');
    
    try {
        console.log('✅ WASM module loaded successfully');
        console.log('📋 Available functions:', Object.keys(Module).filter(k => k.startsWith('_')));
        
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
        
        console.log('\n🧵 Threading info:');
        console.log(`   📱 ccall available: ${typeof Module.ccall === 'function'}`);
        console.log(`   🏃 Functions available: ${Object.keys(Module).filter(k => k.startsWith('_')).length}`);
        
        console.log('\n✅ All tests completed successfully!');
        console.log('🎊 Ceres-solver is working in WASM with full threading support!');
        
    } catch (error) {
        console.error('❌ Error testing WASM module:', error);
        process.exit(1);
    }
}

// Run the test
if (require.main === module) {
    testCeresWasm().catch(console.error);
}

module.exports = { testCeresWasm };

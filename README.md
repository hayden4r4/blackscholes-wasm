# blackscholes  
  
This library provides an simple, lightweight, and efficient (though not heavily optimized) implementation of the Black-Scholes-Merton model for pricing European options.  
  
Includes all first, second, and third order Greeks.  

Implements:  
  
- calc_iv() in the ImpliedVolatility trait which uses [Modified Corrado-Miller by Piotr P√luciennik (2007)](https://sin.put.poznan.pl/files/download/37938) for the initial volatility guess and the Newton Raphson algorithm to solve for the implied volatility.
  
## Usage  
  
View the [docs](https://docs.rs/blackscholes_wasm) for usage and examples.  
  
**Other packages available:**  
Python: [Pypi](https://pypi.org/project/blackscholes-python/)  
Rust: [crates.io](https://crates.io/crates/blackscholes)  

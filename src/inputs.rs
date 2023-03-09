use std::fmt::{Display, Formatter, Result as fmtResult};
use wasm_bindgen::prelude::*;

/// The type of option to be priced (call or put).
#[wasm_bindgen]
#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub enum OptionType {
    Call,
    Put,
}

impl Display for OptionType {
    fn fmt(&self, f: &mut Formatter) -> fmtResult {
        match self {
            OptionType::Call => write!(f, "Call"),
            OptionType::Put => write!(f, "Put"),
        }
    }
}

/// The inputs to the Black-Scholes-Merton model.
#[wasm_bindgen(inspectable)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Inputs {
    /// The type of the option (call or put)
    pub option_type: OptionType,
    /// Stock price
    pub s: f32,
    /// Strike price
    pub k: f32,
    /// Option price
    pub p: Option<f32>,
    /// Risk-free rate
    pub r: f32,
    /// Dividend yield
    pub q: f32,
    /// Time to maturity in years
    pub t: f32,
    /// Volatility
    pub sigma: Option<f32>,
}

/// Methods for calculating the price, greeks, and implied volatility of an option.
#[wasm_bindgen]
impl Inputs {
    /// Creates instance ot the `Inputs` struct.
    /// # Arguments
    /// * `option_type` - The type of option to be priced.
    /// * `s` - The current price of the underlying asset.
    /// * `k` - The strike price of the option.
    /// * `p` - The dividend yield of the underlying asset.
    /// * `r` - The risk-free interest rate.
    /// * `q` - The dividend yield of the underlying asset.
    /// * `t` - The time to maturity of the option in years.
    /// * `sigma` - The volatility of the underlying asset.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// ```
    /// # Returns
    /// An instance of the `Inputs` struct.
    #[wasm_bindgen(constructor)]
    pub fn new(
        option_type: OptionType,
        s: f32,
        k: f32,
        p: Option<f32>,
        r: f32,
        q: f32,
        t: f32,
        sigma: Option<f32>,
    ) -> Self {
        Self {
            option_type,
            s,
            k,
            p,
            r,
            q,
            t,
            sigma,
        }
    }
}

impl Display for Inputs {
    fn fmt(&self, f: &mut Formatter) -> fmtResult {
        writeln!(f, "Option type: {}", self.option_type)?;
        writeln!(f, "Stock price: {:.2}", self.s)?;
        writeln!(f, "Strike price: {:.2}", self.k)?;
        match self.p {
            Some(p) => writeln!(f, "Option price: {:.2}", p)?,
            None => writeln!(f, "Option price: None")?,
        }
        writeln!(f, "Risk-free rate: {:.4}", self.r)?;
        writeln!(f, "Dividend yield: {:.4}", self.q)?;
        writeln!(f, "Time to maturity: {:.4}", self.t)?;
        match self.sigma {
            Some(sigma) => writeln!(f, "Volatility: {:.4}", sigma)?,
            None => writeln!(f, "Volatility: None")?,
        }
        Ok(())
    }
}

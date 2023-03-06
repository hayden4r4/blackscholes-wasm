use crate::{common::*, constants::*, Inputs, OptionType};
use num_traits::Float;
use serde::Serialize;
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
impl Inputs {
    /// Calculates the delta of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the delta of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let delta = inputs.calc_delta().unwrap();
    /// ```
    pub fn calc_delta(&self) -> Result<f32, String> {
        let (nd1, _): (f32, f32) = calc_nd1nd2(&self)?;
        let delta: f32 = match self.option_type {
            OptionType::Call => nd1 * E.powf(-self.q * self.t),
            OptionType::Put => -nd1 * E.powf(-self.q * self.t),
        };
        Ok(delta)
    }

    /// Calculates the gamma of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the gamma of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let gamma = inputs.calc_gamma().unwrap();
    /// ```
    pub fn calc_gamma(&self) -> Result<f32, String> {
        let sigma = self
            .sigma
            .ok_or("Expected Some(f32) for self.sigma, received None")?;

        let nprimed1: f32 = calc_nprimed1(&self)?;
        let gamma: f32 = E.powf(-self.q * self.t) * nprimed1 / (self.s * sigma * self.t.sqrt());
        Ok(gamma)
    }

    /// Calculates the theta of the option.
    /// Uses 365.25 days in a year for calculations.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of theta per day (not per year).
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let theta = inputs.calc_theta().unwrap();
    /// ```
    pub fn calc_theta(&self) -> Result<f32, String> {
        let sigma = self
            .sigma
            .ok_or("Expected Some(f32) for self.sigma, received None")?;

        let nprimed1: f32 = calc_nprimed1(&self)?;
        let (nd1, nd2): (f32, f32) = calc_nd1nd2(&self)?;

        // Calculation uses 365.25 for f32: Time of days per year.
        let theta: f32 = match self.option_type {
            OptionType::Call => {
                (-(self.s * sigma * E.powf(-self.q * self.t) * nprimed1 / (2.0 * self.t.sqrt()))
                    - self.r * self.k * E.powf(-self.r * self.t) * nd2
                    + self.q * self.s * E.powf(-self.q * self.t) * nd1)
                    / DAYS_PER_YEAR
            }
            OptionType::Put => {
                (-(self.s * sigma * E.powf(-self.q * self.t) * nprimed1 / (2.0 * self.t.sqrt()))
                    + self.r * self.k * E.powf(-self.r * self.t) * nd2
                    - self.q * self.s * E.powf(-self.q * self.t) * nd1)
                    / DAYS_PER_YEAR
            }
        };
        Ok(theta)
    }

    /// Calculates the vega of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the vega of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let vega = inputs.calc_vega().unwrap();
    /// ```
    pub fn calc_vega(&self) -> Result<f32, String> {
        let nprimed1: f32 = calc_nprimed1(&self)?;
        let vega: f32 = 0.01 * self.s * E.powf(-self.q * self.t) * self.t.sqrt() * nprimed1;
        Ok(vega)
    }

    /// Calculates the rho of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the rho of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let rho = inputs.calc_rho().unwrap();
    /// ```
    pub fn calc_rho(&self) -> Result<f32, String> {
        let (_, nd2): (f32, f32) = calc_nd1nd2(&self)?;
        let rho: f32 = match &self.option_type {
            OptionType::Call => 1.0 / 100.0 * self.k * self.t * E.powf(-self.r * self.t) * nd2,
            OptionType::Put => -1.0 / 100.0 * self.k * self.t * E.powf(-self.r * self.t) * nd2,
        };
        Ok(rho)
    }

    // The formulas for the greeks below are from the wikipedia page for the Black-Scholes greeks
    // https://en.wikipedia.org/wiki/Greeks_(finance)#Black.E2.80.93Scholes_Greeks
    // Some sources I reviewed contain variations of these formulas and/or varying values, therefore the
    // values returned by this library may not match other libraries or sources.
    // These functions have not been throughouly tested.

    /// Calculates the epsilon of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the epsilon of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let epsilon = inputs.calc_epsilon().unwrap();
    /// ```
    pub fn calc_epsilon(&self) -> Result<f32, String> {
        let (nd1, _) = calc_nd1nd2(&self)?;
        let e_negqt = E.powf(-self.q * self.t);
        let epsilon: f32 = match &self.option_type {
            OptionType::Call => -self.s * self.t * e_negqt * nd1,
            OptionType::Put => self.s * self.t * e_negqt * nd1,
        };
        Ok(epsilon)
    }

    /// Calculates the lambda of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the lambda of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let lambda = inputs.calc_lambda().unwrap();
    /// ```
    pub fn calc_lambda(&self) -> Result<f32, String> {
        let delta = self.calc_delta()?;
        Ok(delta * self.s / self.calc_price()?)
    }

    /// Calculates the vanna of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the vanna of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let vanna = inputs.calc_vanna().unwrap();
    /// ```
    pub fn calc_vanna(&self) -> Result<f32, String> {
        let sigma = self
            .sigma
            .ok_or("Expected Some(f32) for self.sigma, received None")?;

        let nprimed1 = calc_nprimed1(&self)?;
        let (_, d2) = calc_d1d2(&self)?;
        let vanna: f32 = d2 * E.powf(-self.q * self.t) * nprimed1 * -0.01 / sigma;
        Ok(vanna)
    }

    // /// Calculates the charm of the option.
    // /// # Requires
    // /// s, k, r, q, t, sigma
    // /// # Returns
    // /// f32 of the charm of the option.
    // /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let charm = inputs.calc_charm().unwrap();
    /// ```
    pub fn calc_charm(&self) -> Result<f32, String> {
        let sigma = self
            .sigma
            .ok_or("Expected Some(f32) for self.sigma, received None")?;
        let nprimed1 = calc_nprimed1(&self)?;
        let (nd1, _) = calc_nd1nd2(&self)?;
        let (_, d2) = calc_d1d2(&self)?;
        let e_negqt = E.powf(-self.q * self.t);

        let charm = match &self.option_type {
            OptionType::Call => {
                self.q * e_negqt * nd1
                    - e_negqt
                        * nprimed1
                        * (2.0 * (self.r - self.q) * self.t - d2 * sigma * self.t.sqrt())
                        / (2.0 * self.t * sigma * self.t.sqrt())
            }
            OptionType::Put => {
                -self.q * e_negqt * nd1
                    - e_negqt
                        * nprimed1
                        * (2.0 * (self.r - self.q) * self.t - d2 * sigma * self.t.sqrt())
                        / (2.0 * self.t * sigma * self.t.sqrt())
            }
        };
        Ok(charm)
    }

    /// Calculates the veta of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the veta of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let veta = inputs.calc_veta().unwrap();
    /// ```
    pub fn calc_veta(&self) -> Result<f32, String> {
        let sigma = self
            .sigma
            .ok_or("Expected Some(f32) for self.sigma, received None")?;
        let nprimed1 = calc_nprimed1(&self)?;
        let (d1, d2) = calc_d1d2(&self)?;
        let e_negqt = E.powf(-self.q * self.t);

        let veta = -self.s
            * e_negqt
            * nprimed1
            * self.t.sqrt()
            * (self.q + ((self.r - self.q) * d1) / (sigma * self.t.sqrt())
                - ((1.0 + d1 * d2) / (2.0 * self.t)));
        Ok(veta)
    }

    /// Calculates the vomma of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the vomma of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let vomma = inputs.calc_vomma().unwrap();
    /// ```
    pub fn calc_vomma(&self) -> Result<f32, String> {
        let sigma = self
            .sigma
            .ok_or("Expected Some(f32) for self.sigma, received None")?;
        let (d1, d2) = calc_d1d2(&self)?;

        let vomma = Inputs::calc_vega(&self)? * ((d1 * d2) / sigma);
        Ok(vomma)
    }

    /// Calculates the speed of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the speed of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let speed = inputs.calc_speed().unwrap();
    /// ```
    pub fn calc_speed(&self) -> Result<f32, String> {
        let sigma = self
            .sigma
            .ok_or("Expected Some(f32) for self.sigma, received None")?;
        let (d1, _) = calc_d1d2(&self)?;
        let gamma = Inputs::calc_gamma(&self)?;

        let speed = -gamma / self.s * (d1 / (sigma * self.t.sqrt()) + 1.0);
        Ok(speed)
    }

    /// Calculates the zomma of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the zomma of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let zomma = inputs.calc_zomma().unwrap();
    /// ```
    pub fn calc_zomma(&self) -> Result<f32, String> {
        let sigma = self
            .sigma
            .ok_or("Expected Some(f32) for self.sigma, received None")?;
        let (d1, d2) = calc_d1d2(&self)?;
        let gamma = Inputs::calc_gamma(&self)?;

        let zomma = gamma * ((d1 * d2 - 1.0) / sigma);
        Ok(zomma)
    }

    /// Calculates the color of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the color of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let color = inputs.calc_color().unwrap();
    /// ```
    pub fn calc_color(&self) -> Result<f32, String> {
        let sigma = self
            .sigma
            .ok_or("Expected Some(f32) for self.sigma, received None")?;
        let (d1, d2) = calc_d1d2(&self)?;
        let nprimed1 = calc_nprimed1(&self)?;
        let e_negqt = E.powf(-self.q * self.t);

        let color = -e_negqt
            * (nprimed1 / (2.0 * self.s * self.t * sigma * self.t.sqrt()))
            * (2.0 * self.q * self.t
                + 1.0
                + (2.0 * (self.r - self.q) * self.t - d2 * sigma * self.t.sqrt())
                    / (sigma * self.t.sqrt())
                    * d1);
        Ok(color)
    }

    /// Calculates the ultima of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the ultima of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let ultima = inputs.calc_ultima().unwrap();
    /// ```
    pub fn calc_ultima(&self) -> Result<f32, String> {
        let sigma = self
            .sigma
            .ok_or("Expected Some(f32) for self.sigma, received None")?;
        let (d1, d2) = calc_d1d2(&self)?;
        let vega = Inputs::calc_vega(&self)?;

        let ultima =
            -vega / sigma.powf(2.0) * (d1 * d2 * (1.0 - d1 * d2) + d1.powf(2.0) + d2.powf(2.0));
        Ok(ultima)
    }

    /// Calculates the dual delta of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the dual delta of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let dual_delta = inputs.calc_dual_delta().unwrap();
    /// ```
    pub fn calc_dual_delta(&self) -> Result<f32, String> {
        let (_, nd2) = calc_nd1nd2(&self)?;
        let e_negqt = E.powf(-self.q * self.t);

        let dual_delta = match self.option_type {
            OptionType::Call => -e_negqt * nd2,
            OptionType::Put => e_negqt * nd2,
        };
        Ok(dual_delta)
    }

    /// Calculates the dual gamma of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the dual gamma of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let dual_gamma = inputs.calc_dual_gamma().unwrap();
    /// ```
    pub fn calc_dual_gamma(&self) -> Result<f32, String> {
        let sigma = self
            .sigma
            .ok_or("Expected Some(f32) for self.sigma, received None")?;
        let nprimed2 = calc_nprimed2(&self)?;
        let e_negqt = E.powf(-self.q * self.t);

        let dual_gamma = e_negqt * (nprimed2 / (self.k * sigma * self.t.sqrt()));
        Ok(dual_gamma)
    }

    /// Calculates all Greeks of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// HashMap of type <String, f32> of all Greeks of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let greeks = inputs.calc_all_greeks().unwrap();
    /// ```
    pub fn calc_all_greeks(&self) -> Result<String, String> {
        #[derive(Serialize)]
        struct Greeks<T>
        where
            T: Float,
        {
            delta: T,
            gamma: T,
            theta: T,
            vega: T,
            rho: T,
            epsilon: T,
            lambda: T,
            vanna: T,
            charm: T,
            veta: T,
            vomma: T,
            speed: T,
            zomma: T,
            color: T,
            ultima: T,
            dual_delta: T,
            dual_gamma: T,
        }
        let greeks: Greeks<f32> = Greeks {
            delta: self.calc_delta()?,
            gamma: self.calc_gamma()?,
            theta: self.calc_theta()?,
            vega: self.calc_vega()?,
            rho: self.calc_rho()?,
            epsilon: self.calc_epsilon()?,
            lambda: self.calc_lambda()?,
            vanna: self.calc_vanna()?,
            charm: self.calc_charm()?,
            veta: self.calc_veta()?,
            vomma: self.calc_vomma()?,
            speed: self.calc_speed()?,
            zomma: self.calc_zomma()?,
            color: self.calc_color()?,
            ultima: self.calc_ultima()?,
            dual_delta: self.calc_dual_delta()?,
            dual_gamma: self.calc_dual_gamma()?,
        };
        Ok(serde_json::to_string(&greeks)
            .map_err(|e| format!("Failed to serialize Greeks, error: {e}"))?)
    }
}

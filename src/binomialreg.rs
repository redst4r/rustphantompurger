// use std::iter::zip;
use bustools::utils::{self, get_progressbar};
use itertools::izip;
use statrs::distribution::Binomial;
use statrs::distribution::Discrete;
use crate::utils::logspace;

/// Estimating 1-SIHR
/// 
/// # Parameters:
/// * z: the number of non-chimeric molecules at aplification r
/// * mr: is the total number of mulecules at amplification r
/// * r: amplification factor
///
/// Note: this does inference over 1-SIHR, i.e. a very small number, so that logspace makes sense
///  
pub fn phantom_binomial_regression(z: &[usize], mr: &[usize], r:&[usize], s: usize) -> (f64, Vec<f64>, Vec<f64>){

    assert_eq!(z.len() , r.len());

    let s_f64 = s as f64;
    let n = 10000;
    // let n = 100;
    // let mut prange: Vec<f64> = linspace(1e-7, 1.0, n);
    let mut prange: Vec<f64> = logspace(-7.0, 0.0, n);
    // remove the 0
    // prange.remove(0);
    // remove the 1
    prange.pop();

    // add 1: not possible atm, as the likelihood will eval to nan
    // prange.push(1.0);

    let mut loglike_range: Vec<f64> = Vec::with_capacity(n);
    let bar = get_progressbar(n as u64);
    for _p in prange.iter(){

        // each (z,mr,r) tuple corresponds to a z = Binomial(N=mr , p^r)
        let mut loglike = 0.0;
        for (zi, mri, ri) in izip!(z, mr, r){
            loglike += loglike_fn(*_p, *zi, *mri, *ri, s_f64);
        }
        bar.inc(1);
        // println!("{}, {}", p, loglike);
        loglike_range.push(loglike);
    }

    let (ix_max, _loglike_max) = utils::argsort::argmax_float(&loglike_range);
    let pmax = prange[ix_max];
    (pmax, prange, loglike_range)
}


fn loglike_fn(p: f64, zi: usize, mri: usize, ri: usize, s_f64: f64) -> f64{
    let ri_f64 = ri as f64;

    // the SIHR (1-p) isnt directly going into the binomial distribution
    // but a transformed version -> p_binomial
    let correction_factor = (s_f64-1.0)* (p/(s_f64-1.0)).powf(ri_f64);

    let mut p_binomial = (1.0-p).powf(ri_f64) + correction_factor;
    if p_binomial > 1.0{
        p_binomial= 1.0;
    }
    let b_rv = Binomial::new(p_binomial, mri as u64).unwrap_or_else(|_| panic!("p {}", p_binomial));
    let loglike_inc = b_rv.ln_pmf(zi as u64);

            // let loglike_inc = binomial_loglike(*zi, *mri, p_binomial);
            // println!("{},{}", loglike_inc2, loglike_inc);
            // println!("p_binomial {}, ri {}, zi {}, mri {}, logp {}", p_binomial, ri, zi, mri, loglike_inc);
    loglike_inc
}

#[cfg(test)]
mod test{
}
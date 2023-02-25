// use std::iter::zip;
use rustbustools::utils;


fn logfactorial(x: usize) -> f64{
    (1..(x+1)).map(|q|(q as f64).ln()).sum()
}

fn log_binomial_coeff(n: usize, k: usize) -> f64{
    logfactorial(n) - logfactorial(k) - logfactorial(n-k)
}

fn binomial_loglike(x: usize, n: usize, p: f64) -> f64{
    // for a single datapoint
    let logcoeff = log_binomial_coeff(n, x);
    let loglike = (x as f64) * p.ln() + ((n-x) as f64) * (1.0-p).ln() + logcoeff;
    loglike
}

fn linspace(min: f64, max: f64, n: usize) -> Vec<f64>{
    assert!(n>=2);    
    let delta = (max-min)/((n-1) as f64);
    let mut x = Vec::with_capacity(n);
    for i in 0..n{
        let y = min + (i as f64)*delta;
        x.push(y);
    }
    x
}

fn logspace(logmin: f64, logmax:f64, n:usize)-> Vec<f64> {
    let logx = linspace(logmin, logmax, n);
    let x = logx.into_iter().map(|y| 10_f64.powf(y)).collect();
    x
}


#[test]
fn test_linspace(){
    let min = 1.0;
    let max = 10.0;
    let n = 10;
    let x = linspace(min, max, n);
    // println!("{:?}",x);
    assert_eq!(x.len(), n);
    assert_eq!(*x.first().unwrap(), min);
    assert_eq!(*x.last().unwrap(), max);
}

#[test]
fn test_logspace(){
    let logmin = 0.0;  
    let logmax = 1.0;
    let n = 10;
    let x = logspace(logmin, logmax, n);
    // println!("{:?}",x);
    assert_eq!(x.len(), n);
    assert_eq!(*x.first().unwrap(), 1.0);
    assert_eq!(*x.last().unwrap(), 10.0);
}

#[test]
fn test_logspace2(){
    let logmin = -6.0;  
    let logmax = 0.0;
    let n = 7;
    let x = logspace(logmin, logmax, n);
    println!("{:?}",x);
    assert_eq!(x.len(), n);
    assert_eq!(*x.first().unwrap(), 1e-6);
    assert_eq!(*x.last().unwrap(), 1.0);
}



use itertools::izip;

pub fn phantom_binomial_regression(z: &[usize], mr: &[usize], r:&[usize], s: usize) -> (f64, Vec<f64>, Vec<f64>){
    // z is the number of non-chimeric molecules at aplification r
    // mr is the total number of mulecules at amplification r
    // amplification r
    assert_eq!(z.len() , r.len());

    let s_f64 = s as f64;

    let n = 100000;
    // let n = 100;
    // let mut prange: Vec<f64> = linspace(1e-7, 1.0, n);
    let mut prange: Vec<f64> = logspace(-7.0, 0.0, n);
    // remove the 0
    // prange.remove(0);
    // remove the 1
    prange.pop();

    // add 1: not possible atm, as the likelihood will eval to nan
    // prange.push(1.0);

    let mut loglike_range: Vec<f64> = Vec::new();
    for p in prange.iter(){

        // each (z,mr,r) tuple corresponds to a z = Binomial(N=mr , p^r)
        let mut loglike = 0.0;
        for (zi, mri, ri) in izip!(z, mr, r){
            let ri_f64 = *ri as f64;

            let correction_factor = (s_f64-1.0)* ((1.0-p)/(s_f64-1.0)).powf(ri_f64);
            let p_binomial = p.powf(ri_f64) + correction_factor;

            loglike += binomial_loglike(*zi, *mri, p_binomial);
        }
        loglike_range.push(loglike);
    }

    // println!("{:?}", loglike_range);
    let (ix_max, _loglike_max) = utils::argsort::argmax_float(&loglike_range);
    let pmax = prange[ix_max];
    (pmax, prange, loglike_range)
}

#[cfg(test)]
mod test{
    use statrs::assert_almost_eq;

    use super::{logfactorial, log_binomial_coeff, binomial_loglike};
    #[test]
    fn test_logfac(){

        let tolerance = 0.000000001;
        assert_almost_eq!(logfactorial(0), 0.0, tolerance);
        assert_almost_eq!(logfactorial(1), 0.0, tolerance);
        assert_almost_eq!(logfactorial(2), 2_f64.ln(), tolerance);
        assert_almost_eq!(logfactorial(3), 6_f64.ln(), tolerance);
    }

    #[test]
    fn test_log_binomial_coeff(){
        let tolerance = 0.000000001;
        
        assert_almost_eq!(log_binomial_coeff(5, 1), 5_f64.ln(), tolerance);
        assert_almost_eq!(log_binomial_coeff(5, 5), 1_f64.ln(), tolerance);
    }

    #[test]
    fn test_binomial_loglike(){
        assert_eq!(
            binomial_loglike(1, 1, 0.5 ),
            0.5_f64.ln()
        )
    }
}
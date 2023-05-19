use std::collections::HashMap;

pub fn logsumexp(x: &[f64]) -> f64 {
    // logsumexp trick

    // getting the max, stupid f64, cant do .iter().max()
    let c = x.iter().reduce(|a, b| if a > b { a } else { b }).unwrap();
    let mut exp_vec: Vec<f64> = Vec::with_capacity(x.len()); // storing exp(x-c)
    for el in x.iter() {
        exp_vec.push((el - c).exp());
    }
    let z: f64 = c + exp_vec.iter().sum::<f64>().ln();
    z
}

// Takes a Hashmap K->V and transforms all values, leaving the keys as they are
// note that this consumes the hashmap!
pub fn valmap<K, V, V2, F>(fun: F, the_map: HashMap<K, V>) -> HashMap<K, V2>
where
    F: Fn(V) -> V2,
    K: Eq + std::hash::Hash,
{
    let r: HashMap<K, V2> = the_map.into_iter().map(|(k, v)| (k, fun(v))).collect();
    r
}

pub fn valmap_ref<K, V, V2, F>(fun: F, the_map: &HashMap<K, V>) -> HashMap<K, V2>
where
    F: Fn(&V) -> V2,
    K: Eq + std::hash::Hash + Clone,
{
    let r: HashMap<K, V2> = the_map.iter().map(|(k, v)| (k.clone(), fun(v))).collect();
    r
}

#[test]
fn test_valmap() {
    let h: HashMap<&str, usize> = vec![("A", 1), ("B", 2)].into_iter().collect();

    // basics
    let r = valmap(|x| x.to_string(), h);

    assert_eq!(r.get("A").unwrap(), "1");
    assert_eq!(r.get("B").unwrap(), "2");

    // capturing some context
    let h: HashMap<&str, usize> = vec![("A", 1), ("B", 2)].into_iter().collect();
    let inc = 10_usize;
    let r = valmap(|x| x + inc, h);

    assert_eq!(*r.get("A").unwrap(), 11);
    assert_eq!(*r.get("B").unwrap(), 12);
}

// fn logfactorial(x: usize) -> f64{
//     (1..(x+1)).map(|q|(q as f64).ln()).sum()
// }

// fn log_binomial_coeff(n: usize, k: usize) -> f64{
//     // this function is RIDICULOUSLY SLOW!!!
//     logfactorial(n) - logfactorial(k) - logfactorial(n-k)
// }

// pub fn binomial_loglike(x: usize, n: usize, p: f64) -> f64{
//     // for a single datapoint

//     //edge cases
//     if p==1.0{
//         if x==n{
//             return 0.0 // log(p=1)
//         }
//         else{
//             return 0_f64.ln()   //zero prob
//         }
//     }
//     if p==0.0{
//         if x==0{
//             return 0.0  // log(p=1)
//         }
//         else{
//             return 0_f64.ln()   //zero prob
//         }
//     }

//     let logcoeff = log_binomial_coeff(n, x);
//     (x as f64) * p.ln() + ((n-x) as f64) * (1.0-p).ln() + logcoeff

// }

pub fn linspace(min: f64, max: f64, n: usize) -> Vec<f64> {
    assert!(n >= 2);
    let delta = (max - min) / ((n - 1) as f64);
    let mut x = Vec::with_capacity(n);
    for i in 0..n {
        let y = min + (i as f64) * delta;
        x.push(y);
    }
    x
}

pub fn logspace(logmin: f64, logmax: f64, n: usize) -> Vec<f64> {
    let logx = linspace(logmin, logmax, n);
    logx.into_iter().map(|y| 10_f64.powf(y)).collect()
}

#[cfg(test)]
mod tests {

    use crate::utils::{linspace, logspace};

    // use crate::utils::{logfactorial, log_binomial_coeff, binomial_loglike, linspace, logspace, logsumexp};
    // use statrs::assert_almost_eq;

    // #[test]
    // fn test_losumexp(){
    //     assert_eq!(logsumexp(&vec![0.0]), 0.0);
    //     assert_eq!(logsumexp(&vec![10.0]), 10.0);
    //     assert_eq!(logsumexp(&vec![0.0, 0.0, 0.0]), 3_f64.ln());
    //     assert_eq!(logsumexp(&vec![1.0, 1.0, 1.0]), 3_f64.ln() + 1.0);  // log[3 e]
    // }

    // #[test]
    // fn test_logfac(){
    //     let tolerance = 0.000000001;
    //     assert_almost_eq!(logfactorial(0), 0.0, tolerance);
    //     assert_almost_eq!(logfactorial(1), 0.0, tolerance);
    //     assert_almost_eq!(logfactorial(2), 2_f64.ln(), tolerance);
    //     assert_almost_eq!(logfactorial(3), 6_f64.ln(), tolerance);
    // }

    // #[test]
    // fn test_log_binomial_coeff(){
    //     let tolerance = 0.000000001;
    //     assert_almost_eq!(log_binomial_coeff(5, 1), 5_f64.ln(), tolerance);
    //     assert_almost_eq!(log_binomial_coeff(5, 5), 1_f64.ln(), tolerance);
    // }

    // #[test]
    // fn test_binomial_loglike(){
    //     assert_eq!(
    //         binomial_loglike(1, 1, 0.5 ),
    //         0.5_f64.ln()
    //     );

    //     // edge cases
    //     assert_eq!(
    //         binomial_loglike(1, 1, 1.0 ),
    //         1_f64.ln()
    //     );
    //     assert_eq!(
    //         binomial_loglike(0, 1, 0.0 ),
    //         1_f64.ln()
    //     );
    // }
    #[test]
    fn test_linspace() {
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
    fn test_logspace() {
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
    fn test_logspace2() {
        let logmin = -6.0;
        let logmax = 0.0;
        let n = 7;
        let x = logspace(logmin, logmax, n);
        println!("{:?}", x);
        assert_eq!(x.len(), n);
        assert_eq!(*x.first().unwrap(), 1e-6);
        assert_eq!(*x.last().unwrap(), 1.0);
    }
}

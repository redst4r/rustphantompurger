use std::{collections::{HashMap}, hash::Hash};
use crate::{utils::{logsumexp}, phantompurger::{FingerprintHistogram, Fingerprint, groupby_gene_across_samples, make_fingerprint_simple, make_fingerprint_histogram}};
use itertools::{izip};
use rustbustools::{bus_multi::{CellUmiIteratorMulti}, io::{BusIteratorBuffered}, iterators::{CbUmiGroupIterator}};
use rustbustools::io::{BusRecord, BusFolder, BusWriter};
use rustbustools::utils::{get_progressbar, argsort::argmax_float};



#[derive(Eq, PartialEq, Hash, Ord, PartialOrd, Clone, Debug, Copy)]
pub struct AmpFactor(u32);

pub struct PhantomPosterior{
    order: Vec<String>,
    pi_norm: HashMap<(AmpFactor, String), f64>,
    p_no_hop: f64,
    vr_norm: HashMap<(AmpFactor, String), f64>
}


// fn record_dict_to_fingerprints(
//     order: &[String], 
//     record_dict: &HashMap<String, Vec<BusRecord>>, 
//     ecmapper_dict:&HashMap<String, &Ec2GeneMapper>) -> Vec<Vec<u32>>{
//     // turn a record dict (as emitted by the multi_iterators)
//     // into a list of fingerprints
//     // 
//     // usually only on fingerprint comes back (if a CB/UMI maps to the same gene across the samples)
//     // but sometimes it gets split

//     let mut fps: Vec<Vec<u32>> = Vec::new();
//     let grouped_record_dicts = groupby_gene_across_samples(&record_dict, &ecmapper_dict);
//     for rd in grouped_record_dicts{


//         let fp_hash = make_fingerprint_simple(&rd);

//         // turn the hashmap into a vector, sorted acc to order
//         let fp:  Vec<_>  = order.iter()
//             .map(|s| fp_hash.get(s).unwrap_or(&0)).cloned().collect();
//         fps.push(fp);
//     }
//     fps
// }

impl PhantomPosterior{

    pub fn new(fph: &FingerprintHistogram) -> Self{

        let order: Vec<String> = fph.order.to_vec();
        let n_samples = order.len() as f64;

        println!("Estimate SIHR");
        let p = 1.0-fph.estimate_sihr(); // the prop of not hopping
        println!("1-SIHR {p}");

        /*
        Estimate the quantities needed for the posterior:
        vr, vr_norm, pi_norm
         */

        // v_rs: the molecular complexity profile of samples
        // for each sample/Amp factor, count the number of reads
        let mut vr: HashMap<(AmpFactor, String), u32> = HashMap::new();
        for (fp, freq) in fph.histogram.iter(){
            let r = AmpFactor(fp.iter().sum());
            for (i, n_reads) in fp.iter().enumerate(){
                // since we observed many of those (its a histogram)
                let n_reads_total = (*n_reads) * (*freq as u32);
                let key =(r, order[i].clone()); 
                let v = vr.entry(key).or_insert(0);
                *v+= n_reads_total
            }
        }
        // for fixed r, normalize across samples
        let mut norm_constant: HashMap<AmpFactor, u32> = HashMap::new();
        for ((r, _s), f) in vr.iter(){
            let v = norm_constant.entry(*r).or_insert(0);
            *v+= *f;
        }

        let vr_norm: HashMap<(AmpFactor, String), f64> = vr
            .into_iter()
            .map(|((r, s),v)| 
                ((r, s), (v as f64) / (*norm_constant.get(&r).unwrap() as f64)))
            .collect();

        // map_to_file(&vr_norm, "/tmp/vr_norm.csv");

        // convert to pi_rs, Eq10 in Suppl
        let pi_norm: HashMap<(AmpFactor, String), f64> = vr_norm
            .iter()
            .map(|((r, s), f)|
                ((*r, s.to_owned()), (f * (n_samples-1.0) + (p-1.0))/ (n_samples*p -1.0))
            )
            // setting too small values to 1e-6
            .map(|((r, s), f)|
                ((r, s), if f>0.0 {f} else {1e-6})
            )
            .collect();      

        // due to the clipping to 1e-6, we need to renorm such that pi_r ==1 
        // "In empirical data, when this relationship is violated, we can set pi_rs = 10−6 and renormalize pi_r accordingly."

        let mut normfactor: HashMap<AmpFactor,f64> = HashMap::new();
        for ((r, _s), f) in pi_norm.iter(){
            let v = normfactor.entry(*r).or_insert(0.0);
            *v+=f;
        }

        let pi_renorm: HashMap<(AmpFactor, String), f64> = pi_norm.iter()
            .map(|((r, s), f)| 
                ((*r, s.clone()), f/normfactor.get(r).unwrap())
            ).collect();

        // map_to_file(&pi_renorm, "/tmp/pi_renorm.csv");

        // let posterior_cache = HashMap::new();
        PhantomPosterior{
            order,
            pi_norm: pi_renorm,  // for each fingerprint a vector of posterior probs
            p_no_hop: p,
            vr_norm
        }
    }

    fn posterior_internal(&self, fingerprint: &[u32]) -> Vec<f64>{
        // Page 22 in Suppl
        let n_samples = self.order.len() as f64;

        assert!(fingerprint.len()==self.order.len());

        let r = AmpFactor(fingerprint.iter().sum());

        let linear = false;

        if linear{
            panic!("creates overflows!");
            let mut posterior: Vec<f64> = Vec::with_capacity(self.order.len());
            let base = (n_samples-1.0)/((1.0/ self.p_no_hop) -1.0);

            for (sname, y) in izip![self.order.iter(), fingerprint.iter() ]{

                let pi = *self.pi_norm.get(&(r, sname.to_string())).unwrap_or_else(|| panic!("Unknown {} {}", sname,r.0));
                let post = base.powi(*y as i32) * pi;
                posterior.push(post);
            }
            let norm_constant: f64 = posterior.iter().sum();
            let posterior_normed: Vec<f64> = posterior.iter().map(|v| (v/norm_constant)).collect();
            posterior_normed
        }
        else{ // log posterior, preventing overflows

            let mut logposterior: Vec<f64> = Vec::with_capacity(self.order.len());
            let logbase = (n_samples-1.0).ln() - ((1.0/ self.p_no_hop) -1.0).ln();
            for (sname, y) in izip![self.order.iter(), fingerprint.iter() ]{

                let pi = *self.pi_norm.get(&(r, sname.to_string())).unwrap_or_else(|| panic!("Unknown {} {}", sname,r.0));
                let logpi = pi.ln();
                let logpost = (*y as f64)*logbase + logpi;
                logposterior.push(logpost);
            }
            let norm_constant: f64 = logsumexp(&logposterior);
            let logposterior_normed: Vec<f64> = logposterior.iter().map(|v| (v-norm_constant)).collect();

            let posterior_normed: Vec<f64> = logposterior_normed.iter().map(|x| x.exp()).collect();
            posterior_normed
        }
    }

    pub fn posterior(&self, fingerprint: Fingerprint) -> HashMap<String, f64>{

        let mut fp_numeric: Vec<u32> = Vec::with_capacity(self.order.len());
        for s in self.order.iter(){
            
            let y = fingerprint.get(s).unwrap_or(&0);
            fp_numeric.push(*y);
        }
        let posterior_numeric = self.posterior_internal(&fp_numeric);

        let mut posterior_map: HashMap<String, f64> = HashMap::new();

        for (sname,post) in izip![self.order.iter(), posterior_numeric.iter()]{
            posterior_map.insert(sname.clone(), *post);
        }
        posterior_map
    }

    pub fn filter_busfiles(
        &self, input_busfolders: HashMap<String, BusFolder>, 
        output_busfolders: HashMap<String, String>,
        output_removed: HashMap<String, String>,
        posterior_threshold: f64
    ){
        // for each fingerprint a vector of posterior probs
        let mut posterior_cache: HashMap<Vec<u32>, Vec<f64>> = HashMap::new();  

        // create the EC2gene mappers
        let ecmapper_dict = input_busfolders.iter()
        .map(|(samplename, bfolder)|
            (samplename.clone(), &bfolder.ec2gene)   //#todo remove clone
        ).collect();

        // a list of busfile-names for the iterator
        // let busnames = valmap(|bfolder| bfolder.get_busfile(), input_busfolders);
        // let busnames = valmap_ref(|bfolder| bfolder.get_busfile(), &input_busfolders);

        let busnames: HashMap<String,String> =input_busfolders.iter()
            .map(|(s, bfolder)| (s.clone(), bfolder.get_busfile()))
            .collect();
    
        let mut buswriters: HashMap<String,BusWriter> = output_busfolders.iter()
            .map(|(sname, fname)| 
                (
                    sname.to_owned(), 
                    BusWriter::new(
                        fname, 
                        input_busfolders.get(sname).unwrap().get_bus_header()
                    )
                ) 
            )
            .collect();

        let mut buswriters_removed: HashMap<String,BusWriter> = output_removed.iter()
            .map(|(sname, fname)| 
                (
                    sname.to_owned(), 
                    BusWriter::new(
                        fname, 
                        input_busfolders.get(sname).unwrap().get_bus_header()
                    )
                ) 
            )
            .collect();

        // figure out size of iterators, just for progress bar!
        let mut total = 0;
        // TODO: this doesnt check if the EC overlaps
        for v in busnames.values(){
            println!("determine size of iterator");
            let total_records = BusIteratorBuffered::new(v).groupby_cbumi().count();
            if total< total_records{
                total=total_records
            }
        }
        println!("total records {}", total);
        let bar = get_progressbar(total as u64);

        let multi_iter = CellUmiIteratorMulti::new(&busnames);
        let mut records_resolved = 0;
        let mut records_ambiguous = 0;
        let mut records_total = 0;
        // let mut records_multi_finger = 0;

        let mut records_original: HashMap<String, usize> = HashMap::new();
        let mut records_written: HashMap<String, usize> = HashMap::new();
        let mut records_filtered: HashMap<String, usize> = HashMap::new();

        let mut reads_original: HashMap<String, usize> = HashMap::new();
        let mut reads_written: HashMap<String, usize> = HashMap::new();
        let mut reads_filtered: HashMap<String, usize> = HashMap::new();

        for (_i,(_cb_umi, record_dict)) in multi_iter.enumerate(){
            records_total += 1;

            if _i % 1_000_000== 0{
                bar.inc(1_000_000);
            }

            // counting
            for (sname, rvector) in record_dict.iter(){
                let v = records_original.entry(sname.to_string()).or_insert(0);
                *v+=rvector.len();

                let v = reads_original.entry(sname.to_string()).or_insert(0);
                let nreads: usize = rvector.iter().map(|r|r.COUNT as usize).sum();
                *v += nreads;
            }

            // turn into consistent cb/umi/gene over samples
            let grouped_record_dicts = groupby_gene_across_samples(&record_dict, &ecmapper_dict);
            for rd in grouped_record_dicts{

                let fp_hash = make_fingerprint_simple(&rd);
        
                // turn the hashmap into a vector, sorted acc to order
                let fp:  &Vec<_>  = &self.order.iter()
                    .map(|s| fp_hash.get(s).unwrap_or(&0)).cloned().collect();

                if !posterior_cache.contains_key(fp){
                    let pvec = self.posterior_internal(fp);
                    posterior_cache.insert(fp.clone(), pvec);
                }
                let posterior_vec = posterior_cache.get(fp).unwrap();

                let (ix, pmax) = argmax_float(posterior_vec);
                let sample_max = self.order[ix].clone();
                // if we find an unambiguous assignment, write the molecule to that sample
                // otherwise drop the molecule in all samples
                if pmax > posterior_threshold {
                    let wr = buswriters.get_mut(&sample_max).unwrap();
                    let r = record_dict.get(&sample_max).unwrap();
                    wr.write_records(r);
                    records_resolved+=1;


                    let v = records_written.entry(sample_max.to_string()).or_insert(0);
                    *v+=r.len();
    
                    let v = reads_written.entry(sample_max.to_string()).or_insert(0);
                    let nreads: usize = r.iter().map(|r|r.COUNT as usize).sum();
                    *v += nreads;
                    

                    // write the filtered reads into the "remove" files
                    for s in self.order.iter(){
                        if *s != sample_max{
                            // if that sample has the molecules, its a phantom
                            if let Some(r) = record_dict.get(s){
                                let wr = buswriters_removed.get_mut(s).unwrap();
                                wr.write_records(r); 


                                let v = records_filtered.entry(s.to_string()).or_insert(0);
                                *v+=r.len();
                
                                let v = reads_filtered.entry(s.to_string()).or_insert(0);
                                let nreads: usize = r.iter().map(|r|r.COUNT as usize).sum();
                                *v += nreads;
                            }
                        }
                    }
                }
                else{
                    // couldnt find clear source sample, just write them back untouched
                    records_ambiguous+=1;
                    // for s in self.order.iter(){
                    //     if let Some(r) = record_dict.get(s){
                    //         let wr = buswriters_removed.get_mut(s).unwrap();
                    //         wr.write_records(r); 
                    //     }
                    // }
                }
            }
        }
        println!("{records_resolved}/{records_total} records corrected/written");
        println!("{records_ambiguous}/{records_total} records were ambigous");
        // println!("{records_multi_finger}/{records_total} records were mutli");


        for s in self.order.iter(){
            let ro = records_original.get(s).unwrap_or(&0);
            let rr =records_written.get(s).unwrap_or(&0);
            let rf = records_filtered.get(s).unwrap_or(&0);
            let total = rr+rf;
            println!(
                "{s} records: orig: {ro}, written {rr}, filtered {rf}, summed: {total}"
            )
        }
        println!("=============================");
        for s in self.order.iter(){
            let ro = reads_original.get(s).unwrap_or(&0);
            let rr =reads_written.get(s).unwrap_or(&0);
            let rf = reads_filtered.get(s).unwrap_or(&0);
            let total = rr+rf;
            println!(
                "{s} reads: orig: {ro}, written {rr}, filtered {rf}, summed: {total}"
            )
        }


        // writing the posterior cache into a file for debug

        // let mut fh = File::create("/tmp/posterior.csv").unwrap();
        // let mut header = self.order.join(",");
        // let header2 = self.order.iter().map(|x| format!("P_{}", x)).join(",");

        // header.push(',');
        // header.push_str(&header2);
        // writeln!(fh, "{}", header).unwrap();
        // for (fp, post) in posterior_cache.iter(){

        //     let fp_string = fp.iter().map(|x| x.to_string()).join(",");
        //     let post_string = post.iter().map(|x| x.to_string()).join(",");
        //     writeln!(fh, "{},{}", fp_string, post_string).unwrap();
        // }

    } 
}   

#[cfg(test)]
mod testing{
    use std::collections::{HashMap};
    use rustbustools::{io::{BusRecord, BusFolder, setup_busfile}};
    use crate::{posterior::PhantomPosterior, phantompurger::{create_dummy_ec, make_fingerprint_histogram}};


    // #[test]
    // fn test_posterior(){
    //     let fph = FingerprintHistogram::from_csv("/home/michi/Dropbox/rustphantompurger/IR56_57_phantom.csv");
    //     let posterior = PhantomPosterior::new(&fph);
    //     map_to_file(&posterior.pi_norm, "/tmp/pi_norm.csv");
    //     map_to_file(&posterior.vr_norm, "/tmp/vr_norm.csv");
    // }
    #[test]
    fn test_filter(){

        // exclusive in 1
        let r1 = BusRecord{CB: 0, UMI: 0, EC: 0, COUNT: 100, FLAG: 0};
        // exclusuve in 2
        let s1 = BusRecord{CB: 0, UMI: 21, EC: 1, COUNT: 100, FLAG: 0};
        // exclusuve in 3
        let t1 = BusRecord{CB: 0, UMI: 1, EC: 0, COUNT: 100, FLAG: 0};

        // shared [10,1,1] in 1/2/3
        let r2 = BusRecord{CB: 1, UMI: 0, EC: 0, COUNT: 100, FLAG: 0};
        let s2 = BusRecord{CB: 1, UMI: 0, EC: 2, COUNT: 1, FLAG: 0};
        let t2 = BusRecord{CB: 1, UMI: 0, EC: 2, COUNT: 1, FLAG: 0};

        let (_b1, tdir1) = setup_busfile(&vec![r1,r2]);
        let (_b2, tdir2) = setup_busfile(&vec![s1,s2]);
        let (_b3, tdir3) = setup_busfile(&vec![t1,t2]);


        let es1 = create_dummy_ec();
        let es2 = create_dummy_ec();
        let es3 = create_dummy_ec();
        let es4 = create_dummy_ec();
        let es5 = create_dummy_ec();
        let es6= create_dummy_ec();
        let d1 = tdir1.path().to_str().unwrap().to_string();
        let d2 = tdir2.path().to_str().unwrap().to_string();
        let d3 = tdir3.path().to_str().unwrap().to_string();

        let busfolders = HashMap::from([
            ("s1".to_string(), BusFolder{foldername: d1.clone(), ec2gene: es1}),
            ("s2".to_string(), BusFolder{foldername: d2.clone(), ec2gene: es2}),
            ("s3".to_string(), BusFolder{foldername: d3.clone(), ec2gene: es3})
        ]);
        let fph = make_fingerprint_histogram(busfolders);
        println!("{:?}", fph.histogram);

        let post = PhantomPosterior::new(&fph);

        let filtered_bus = HashMap::from([
            ("s1".to_string(), "/tmp/s1_filtered.bus".to_string()),
            ("s2".to_string(), "/tmp/s2_filtered.bus".to_string()),
            ("s3".to_string(), "/tmp/s3_filtered.bus".to_string())
        ]);
        let removed_bus = HashMap::from([
            ("s1".to_string(), "/tmp/s1_removed.bus".to_string()),
            ("s2".to_string(), "/tmp/s2_removed.bus".to_string()),
            ("s3".to_string(), "/tmp/s3_removed.bus".to_string()),
        ]);        
        let busfolders = HashMap::from([
            ("s1".to_string(), BusFolder{foldername: d1, ec2gene: es4}),
            ("s2".to_string(), BusFolder{foldername: d2, ec2gene: es5}),
            ("s3".to_string(), BusFolder{foldername: d3, ec2gene: es6}),
        ]);
        post.filter_busfiles(
            busfolders, 
            filtered_bus, 
            removed_bus, 
            0.99
        );
        println!("============================");
        rustbustools::inspect::inspect("/tmp/s1_filtered.bus");
        println!("============================");
        rustbustools::inspect::inspect("/tmp/s2_filtered.bus");
        println!("============================");
        rustbustools::inspect::inspect("/tmp/s3_filtered.bus");


    }

    #[test]
    fn test_fullbus(){

        // let fname1 = "/home/michi/bus_testing/bus_output/";
        let fname1 = "/home/michi/bus_testing/bus_output_short";
        let fname2 = "/home/michi/bus_testing/bus_output_short";

        // let fname1 = "/home/michi/mounts/TB4drive/ISB_data/MNGZ01/MS_processed/S1/kallisto/sort_bus/bus_output/";
        // let fname2 = "/home/michi/mounts/TB4drive/ISB_data/MNGZ01/MS_processed/S2/kallisto/sort_bus/bus_output/";
        
        // let fname1 = "/home/michi/mounts/TB4drive/ISB_data/LT_pilot/LT_pilot/kallisto_quant/Fresh1/kallisto/sort_bus/bus_output/";
        // let fname2 = "/home/michi/mounts/TB4drive/ISB_data/LT_pilot/LT_pilot/kallisto_quant/Fresh2/kallisto/sort_bus/bus_output/";

        // "/home/mstrasse/TB4/resources/transcripts_to_genes.txt"
        let t2g = "/home/michi/mounts/TB4drive/kallisto_resources/transcripts_to_genes.txt";

        let busfolders = HashMap::from([
            ("s1".to_string(), BusFolder::new(fname1, t2g)),
            ("s2".to_string(), BusFolder::new(fname2.clone(),t2g)),
        ]);

        let fph = make_fingerprint_histogram(busfolders);
        let post = PhantomPosterior::new(&fph);

        let filtered_bus = HashMap::from([
            ("s1".to_string(), "/tmp/s1_filtered.bus".to_string()),
            ("s2".to_string(), "/tmp/s2_filtered.bus".to_string()),
        ]);
        let removed_bus = HashMap::from([
            ("s1".to_string(), "/tmp/s1_removed.bus".to_string()),
            ("s2".to_string(), "/tmp/s2_removed.bus".to_string()),
        ]);        
        let busfolders = HashMap::from([
            ("s1".to_string(), BusFolder::new(fname1, t2g)),
            ("s2".to_string(), BusFolder::new(fname2.clone(),t2g)),
        ]);
        post.filter_busfiles(
            busfolders, 
            filtered_bus, 
            removed_bus, 
            0.99
        );

        println!("=========Unfiltered===================");
        rustbustools::inspect::inspect(&format!("{fname1}/output.corrected.sort.bus"));
        println!("===========filtered=================");
        rustbustools::inspect::inspect("/tmp/s1_filtered.bus");

        println!("#########################################");

        println!("=========Unfiltered===================");
        rustbustools::inspect::inspect(&format!("{fname2}/output.corrected.sort.bus"));
        println!("===========filtered=================");
        rustbustools::inspect::inspect("/tmp/s2_filtered.bus");

        // let bar = get_progressbar((b1_records+b2_records) as u64);
        // let mut counter = 0;
        // let mut record_counter: HashMap<String, usize> = HashMap::new();
        // let mut cbumi_counter: HashMap<String, usize> = HashMap::new();

    }
}
use std::{collections::HashMap, fs::File, hash::Hash};
use crate::{phantompurger::{groupby_gene_across_samples, make_fingerprint_simple, Fingerprint, FingerprintHistogram}, utils::{ec_mapper_dict_from_busfolders, logsumexp, valmap_ref}};
use itertools::{izip, Itertools};
use bustools::{consistent_genes::Ec2GeneMapper, io::{BusReader, BusWriterPlain}, iterators::CbUmiGroupIterator, merger::MultiIterator, utils::get_spinner};
use bustools::io::BusFolder;
use bustools::utils::argsort::argmax_float;
use serde;
use std::io::Write;

#[derive(Eq, PartialEq, Hash, Ord, PartialOrd, Clone, Debug, Copy, serde::Serialize)]
pub struct AmpFactor(u32);

pub struct PhantomPosterior{
    samplenames: Vec<String>,
    pi_norm: HashMap<(AmpFactor, String), f64>,
    p_no_hop: f64,
    posterior_cache: HashMap<Vec<u32>, Vec<f64>>  
    // vr_norm: HashMap<(AmpFactor, String), f64>
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

    /// Create a PhantomPosterior from a FingerPrint Histogram
    /// 
    /// Estimates the index hopping rate and other key quantities needed
    /// for posterior assignment of bus records
    pub fn new(fph: &FingerprintHistogram) -> Self{

        let order: Vec<String> = fph.samplenames.to_vec();
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
                *v+= n_reads_total;
            }
        }
        // for fixed r, normalize across samples
        let mut norm_constant: HashMap<AmpFactor, u32> = HashMap::new();
        for ((r, _s), f) in &vr{
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

        let posterior_cache = HashMap::new();
        PhantomPosterior{
            samplenames: order,
            pi_norm: pi_renorm,  // for each fingerprint a vector of posterior probs
            p_no_hop: p,
            posterior_cache
            // vr_norm
        }
    }

    /// calculates the posterior probability of a fingerprint,
    /// i.e. the origin of those reads
    fn posterior_internal(&self, fingerprint: &[u32]) -> Vec<f64>{
        // Page 22 in Suppl
        let n_samples = self.samplenames.len() as f64; // `S` in the suppl

        assert!(fingerprint.len()==self.samplenames.len()); 

        let r = AmpFactor(fingerprint.iter().sum()); // `r`

        let mut logposterior: Vec<f64> = Vec::with_capacity(self.samplenames.len());
        // essentially: fingerprint * log(factor) + log(conditional[r])

        let logbase = (n_samples-1.0).ln() - ((1.0/ self.p_no_hop) -1.0).ln();
        for (sname, y) in izip![self.samplenames.iter(), fingerprint.iter() ]{

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

    /// calculates the posterior assignment probability of a record with a given fingerprint across samples
    /// (intuitively, the read would go to the sampel where it has the highest reads)
    pub fn posterior(&self, fingerprint: Fingerprint) -> HashMap<String, f64>{

        let mut fp_numeric: Vec<u32> = Vec::with_capacity(self.samplenames.len());
        for s in &self.samplenames{
            let y = fingerprint.get(s).unwrap_or(&0);
            fp_numeric.push(*y);
        }
        let posterior_numeric = self.posterior_internal(&fp_numeric);

        let mut posterior_map: HashMap<String, f64> = HashMap::new();

        for (sname,post) in izip![self.samplenames.iter(), posterior_numeric.iter()]{
            posterior_map.insert(sname.clone(), *post);
        }
        posterior_map
    }

    #[allow(non_snake_case)]
    fn find_MAP_sample(&mut self, fp: &Vec<u32>, posterior_threshold: f64) -> Option<String>{

        if !self.posterior_cache.contains_key(fp){
            let pvec = self.posterior_internal(fp);
            self.posterior_cache.insert(fp.clone(), pvec);
        }
        let posterior_vec = self.posterior_cache.get(fp).unwrap();

        let (ix, pmax) = argmax_float(posterior_vec);
        let sample_max = self.samplenames[ix].clone();

        // TODO: slight adaptation: 
        // if the fingerprint is only a single sample
        // lower the threshold
        // If the read is only seen in a single sample
        // it doesnt make alot of sense to assign it somewhere else
        let mut threshold = posterior_threshold;

        let n_samples = fp.iter().filter(|&x| *x>0).count();
        if n_samples == 1{
            threshold = 0.95;
        }

        if pmax > threshold {
            Some(sample_max)
        }
        else{
            None
        }
    }

    /// Filter a collection of Bus-quantifications, removing reads/busrecords
    /// that hopped across samples during sequencing
    /// # Parameters
    /// * input_busfolders: HashMap of Samplename -> BusFolder; all the samples to filter
    /// * output_busfolders: where to write the filtered busfiles
    /// * output_removed: where to write the busrecords that got filtered per sample
    /// * output_ambiguous: where to write the busrecords that were ambigous (cant be assigned to a single sample)
    /// * posterior_threshold: Any record that whose posterior assignement probabilty to a single sample exeeds the thresold is written into the output. Otherwise it goes to ambiguois
    pub fn filter_busfiles(
        &mut self, input_busfolders: &HashMap<String, BusFolder>, 
        output_busfolders: &HashMap<String, String>,
        output_removed: &HashMap<String, String>,
        output_ambiguous: &HashMap<String, String>,
        posterior_threshold: f64,
        t2g_file: &str
    ){
        // create the EC2gene mappers
        // silly, cant create one where ECMapper is a reference
        // need to instantiate/own it, then create ref
        let ec_tmp = ec_mapper_dict_from_busfolders(input_busfolders, t2g_file );
        let ecmapper_dict: HashMap<String, &Ec2GeneMapper> = ec_tmp.iter().map(|(name, d)| (name.clone(), d )).collect();


        // a list of busfile-names for the iterator
        // let busnames = valmap(|bfolder| bfolder.get_busfile(), input_busfolders);
        // let busnames = valmap_ref(|bfolder| bfolder.get_busfile(), &input_busfolders);

        let busnames: HashMap<String,String> =input_busfolders.iter()
            .map(|(s, bfolder)| (s.clone(), bfolder.get_busfile()))
            .collect();
    
        let mut buswriters: HashMap<String,BusWriterPlain> = output_busfolders.iter()
            .map(|(sname, fname)| 
                (
                    sname.to_owned(), 
                    BusWriterPlain::new(
                        fname, 
                        input_busfolders.get(sname).unwrap().get_bus_params()
                    )
                ) 
            )
            .collect();

        let mut buswriters_removed: HashMap<String,BusWriterPlain> = output_removed.iter()
            .map(|(sname, fname)| 
                (
                    sname.to_owned(), 
                    BusWriterPlain::new(
                        fname, 
                        input_busfolders.get(sname).unwrap().get_bus_params()
                    )
                ) 
            )
            .collect();

        let mut buswriters_ambigous: HashMap<String,BusWriterPlain> = output_ambiguous.iter()
            .map(|(sname, fname)| 
                (
                    sname.to_owned(), 
                    BusWriterPlain::new(
                        fname, 
                        input_busfolders.get(sname).unwrap().get_bus_params()
                    )
                ) 
            )
            .collect();

        let bar = get_spinner();

        let iterators = valmap_ref(
            |busfile|{
                BusReader::new(busfile).groupby_cbumi()
            },
            &busnames);

        let multi_iter = MultiIterator::new(iterators);
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
            for (sname, rvector) in &record_dict{
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
                let fp:  &Vec<_>  = &self.samplenames.iter()
                    .map(|s| fp_hash.get(s).unwrap_or(&0)).cloned().collect();

                // if we find an unambiguous assignment, write the molecule to that sample
                // otherwise drop the molecule in all samples
                if let Some(sample_max) = self.find_MAP_sample(fp, posterior_threshold){
                    let wr = buswriters.get_mut(&sample_max).unwrap();
                    let r = rd.get(&sample_max).unwrap();
                    wr.write_record(r);
                    records_resolved+=1;

                    let v = records_written.entry(sample_max.to_string()).or_insert(0);
                    *v+=1;
    
                    let v = reads_written.entry(sample_max.to_string()).or_insert(0);
                    *v += r.COUNT as usize;
                    
                    // write the filtered reads into the "remove" files
                    for s in &self.samplenames{
                        if *s != sample_max{
                            // if that sample has the molecules, its a phantom
                            if let Some(r) = rd.get(s){
                                let wr = buswriters_removed.get_mut(s).unwrap();
                                wr.write_record(r); 

                                let v = records_filtered.entry(s.to_string()).or_insert(0);
                                *v+=1;
                
                                let v = reads_filtered.entry(s.to_string()).or_insert(0);
                                *v += r.COUNT as usize;
                            }
                        }
                    }
                }
                else{
                    for (sname, record) in rd{
                        let wr = buswriters_ambigous.get_mut(&sname).unwrap();
                        wr.write_record(&record); 

                    }
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


        for s in &self.samplenames{
            let ro = records_original.get(s).unwrap_or(&0);
            let rr =records_written.get(s).unwrap_or(&0);
            let rf = records_filtered.get(s).unwrap_or(&0);
            let total = rr+rf;
            println!(
                "{s} records: orig: {ro}, written {rr}, filtered {rf}, summed: {total}"
            );
        }
        println!("=============================");
        for s in &self.samplenames{
            let ro = reads_original.get(s).unwrap_or(&0);
            let rr =reads_written.get(s).unwrap_or(&0);
            let rf = reads_filtered.get(s).unwrap_or(&0);
            let total = rr+rf;
            println!(
                "{s} reads: orig: {ro}, written {rr}, filtered {rf}, summed: {total}"
            );
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

fn map_to_file(h: &HashMap<(AmpFactor, String), f64>, outfile: &str) {
    let mut amp: Vec<AmpFactor> = h.keys().map(|(r, _s)| *r).unique().collect();
    let mut samplenames: Vec<String> = h.keys().map(|(_r, s)| s.to_owned()).unique().collect();

    amp.sort();
    samplenames.sort();

    let mut fh = File::create(outfile).unwrap();
    let mut header = "samplename,".to_string();
    let rstring = amp
        .iter()
        .map(|x| x.0.to_string())
        .collect::<Vec<String>>()
        .join(",");
    header.push_str(&rstring);
    header.push_str(",frequency");

    writeln!(fh, "{}", header).unwrap();

    for s in samplenames {
        let pi_vec: Vec<String> = amp
            .iter()
            .map(|r| h.get(&(*r, s.clone())).unwrap().to_string())
            .collect();
        writeln!(fh, "{},{}", s, pi_vec.join(",")).unwrap();
    }
}

#[cfg(test)]
mod testing{
    use insta;
    use super::*;    
    use crate::posterior::PhantomPosterior;

    #[test]
    fn test_posterior_pinorm(){
        let fph = FingerprintHistogram::from_csv("/home/michi/Dropbox/rustphantompurger/IR56_57_phantom.csv");
        let posterior = PhantomPosterior::new(&fph);
        // map_to_file(&posterior.pi_norm, "/tmp/pi_norm.csv");
        // map_to_file(&posterior.vr_norm, "/tmp/vr_norm.csv");

        insta::with_settings!({sort_maps => true}, {
            insta::assert_yaml_snapshot!(
                posterior.pi_norm,
                {".*" => insta::rounded_redaction(6)}
            );
        });  

    }

    #[test]
    fn test_posterior2(){
        let fph = FingerprintHistogram::from_csv("/home/michi/Dropbox/rustphantompurger/IR56_57_phantom.csv");
        let posterior = PhantomPosterior::new(&fph);
        
        // check that the SIHR is still the same
        // deprecated, already testing in phnatompruger.rs
        let sihr = 1.0-posterior.p_no_hop;
        insta::assert_yaml_snapshot!(sihr, @r###"
        ---
        0.0016278754068794754
        "###);

        let fp: HashMap<String, u32> = vec![
            ("DT75-HL60".to_string(), 1_u32),
            ("DT75-day2".to_string(), 0_u32),
            ("DT75-day3".to_string(), 0_u32),
            ("DT76-HL60".to_string(), 0_u32),
            ("DT76-day1".to_string(), 0_u32),
            ("DT76-day2".to_string(), 0_u32),
            ("DT76-day3".to_string(), 0_u32),
        ].into_iter().collect();

        let post = posterior.posterior(fp);
        // map_to_file(&posterior.vr_norm, "/tmp/vr_norm.csv");


        let fp2: HashMap<String, u32> = vec![
            ("DT75-HL60".to_string(), 1_u32),
            ("DT75-day2".to_string(), 0_u32),
            ("DT75-day3".to_string(), 1_u32),
            ("DT76-HL60".to_string(), 0_u32),
            ("DT76-day1".to_string(), 0_u32),
            ("DT76-day2".to_string(), 0_u32),
            ("DT76-day3".to_string(), 0_u32),
        ].into_iter().collect();
        let post2 = posterior.posterior(fp2);

        // need to enforce hashmap sorting! otherwise the resulting Hashmap (Even if same)
        // will differ from the snapshot on disk
        insta::with_settings!({sort_maps => true}, {
            insta::assert_yaml_snapshot!("posterior_test1", post, {".*" => insta::rounded_redaction(6)});
            insta::assert_yaml_snapshot!("posterior_test2", post2, {".*" => insta::rounded_redaction(6)});
        });  
    }

    /*
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

        let r3 = BusRecord{CB: 1, UMI: 0, EC: 0, COUNT: 100, FLAG: 0};
        let s3 = BusRecord{CB: 1, UMI: 0, EC: 0, COUNT: 1, FLAG: 0};
        let t3 = BusRecord{CB: 1, UMI: 0, EC: 1, COUNT: 1, FLAG: 0};

        let (_b1, tdir1) = setup_busfile(&vec![r1,r2, r3]);
        let (_b2, tdir2) = setup_busfile(&vec![s1,s2, s3]);
        let (_b3, tdir3) = setup_busfile(&vec![t1,t2, t3]);


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

        let mut post = PhantomPosterior::new(&fph);

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
        let ambigous_bus = HashMap::from([
            ("s1".to_string(), "/tmp/s1_ambigous.bus".to_string()),
            ("s2".to_string(), "/tmp/s2_ambigous.bus".to_string()),
            ("s3".to_string(), "/tmp/s3_ambigous.bus".to_string()),
        ]); 
        let busfolders = HashMap::from([
            ("s1".to_string(), BusFolder{foldername: d1, ec2gene: es4}),
            ("s2".to_string(), BusFolder{foldername: d2, ec2gene: es5}),
            ("s3".to_string(), BusFolder{foldername: d3, ec2gene: es6}),
        ]);
        post.filter_busfiles(
            &busfolders, 
            &filtered_bus, 
            &removed_bus, 
            &ambigous_bus,
            0.99
        );
        println!("============================");
        // rustbustools_cli::inspect::inspect("/tmp/s1_filtered.bus");
        // println!("============================");
        // rustbustools_cli::inspect::inspect("/tmp/s2_filtered.bus");
        // println!("============================");
        // rustbustools_cli::inspect::inspect("/tmp/s3_filtered.bus");

    }

    #[test]
    fn test_fullbus(){

        // let fname1 = "/home/michi/bus_testing/bus_output/";
        // let fname1 = "/home/michi/bus_testing/bus_output_short";
        // let fname2 = "/home/michi/bus_testing/bus_output_short";

        // let fname1 = "/home/michi/mounts/TB4drive/ISB_data/MNGZ01/MS_processed/S1/kallisto/sort_bus/bus_output/";
        // let fname2 = "/home/michi/mounts/TB4drive/ISB_data/MNGZ01/MS_processed/S2/kallisto/sort_bus/bus_output/";
        
        let fname1 = "/home/michi/mounts/TB4drive/ISB_data/LT_pilot/LT_pilot/kallisto_quant/Fresh1/kallisto/sort_bus/bus_output/";
        let fname2 = "/home/michi/mounts/TB4drive/ISB_data/LT_pilot/LT_pilot/kallisto_quant/Fresh2/kallisto/sort_bus/bus_output/";

        // "/home/mstrasse/TB4/resources/transcripts_to_genes.txt"
        let t2g = "/home/michi/mounts/TB4drive/kallisto_resources/transcripts_to_genes.txt";

        let busfolders = HashMap::from([
            ("s1".to_string(), BusFolder::new(fname1, t2g)),
            ("s2".to_string(), BusFolder::new(fname2.clone(),t2g)),
        ]);

        let fph = make_fingerprint_histogram(busfolders);
        let mut post = PhantomPosterior::new(&fph);

        let filtered_bus = HashMap::from([
            ("s1".to_string(), "/tmp/s1_filtered.bus".to_string()),
            ("s2".to_string(), "/tmp/s2_filtered.bus".to_string()),
        ]);
        let removed_bus = HashMap::from([
            ("s1".to_string(), "/tmp/s1_removed.bus".to_string()),
            ("s2".to_string(), "/tmp/s2_removed.bus".to_string()),
        ]);       
        let amb_bus = HashMap::from([
            ("s1".to_string(), "/tmp/s1_ambigous.bus".to_string()),
            ("s2".to_string(), "/tmp/s2_ambigous.bus".to_string()),
        ]);  
        let busfolders = HashMap::from([
            ("s1".to_string(), BusFolder::new(fname1, t2g)),
            ("s2".to_string(), BusFolder::new(fname2.clone(),t2g)),
        ]);
        post.filter_busfiles(
            &busfolders, 
            &filtered_bus, 
            &removed_bus, 
            &amb_bus,
            0.99
        );

        // println!("=========Unfiltered===================");
        // rustbustools_cli::inspect::inspect(&format!("{fname1}/output.corrected.sort.bus"));
        // println!("===========filtered=================");
        // rustbustools_cli::inspect::inspect("/tmp/s1_filtered.bus");

        // println!("#########################################");

        // println!("=========Unfiltered===================");
        // rustbustools_cli::inspect::inspect(&format!("{fname2}/output.corrected.sort.bus"));
        // println!("===========filtered=================");
        // rustbustools_cli::inspect::inspect("/tmp/s2_filtered.bus");

        // let bar = get_progressbar((b1_records+b2_records) as u64);
        // let mut counter = 0;
        // let mut record_counter: HashMap<String, usize> = HashMap::new();
        // let mut cbumi_counter: HashMap<String, usize> = HashMap::new();

    }
    */
}

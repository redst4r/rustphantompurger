use std::{collections::{HashMap, HashSet}, fs::File, io::{Write, BufReader, BufRead}, time::Instant, hash::Hash};
use crate::{binomialreg::phantom_binomial_regression, utils::{logsumexp, valmap, valmap_ref}};
use itertools::{izip, Itertools};
use rustbustools::{bus_multi::{CellUmiIteratorMulti, CellIteratorMulti}, io::BusIteratorBuffered, iterators::{CbUmiGroupIterator, CellGroupIterator}};
use rustbustools::io::{BusRecord, BusFolder, BusWriter};
use rustbustools::utils::{get_progressbar, argsort::argmax_float};
use rustbustools::consistent_genes::Ec2GeneMapper;
use rustbustools::disjoint::DisjointSubsets;


#[derive(Eq, PartialEq, Hash, Ord, PartialOrd, Clone, Debug,Copy)]
pub struct CB(u64);

pub fn detect_cell_overlap(busfolders: HashMap<String, String>, outfile: &str) {
    // for each CB calcualte the number of UMIs per experiment/busfile
    // if extensive swapping has occured, we'll see shared CBs across experiments
    // with correlated #umis

    // figure out size of iterators, just for progress bar!
    let mut total = 0;
    // TODO: this doesnt check if the EC overlaps
    for v in busfolders.values(){
        println!("determine size of iterator");
        let total_records = BusIteratorBuffered::new(v).groupby_cb().count();
        if total< total_records{
            total=total_records
        }
    }
    println!("total records {}", total);

    let samplenames: Vec<String> = busfolders.keys().cloned().collect();
    let multi_iter = CellIteratorMulti::new(&busfolders);
    let mut result: HashMap<CB, Vec<usize>> = HashMap::new();

    let bar = get_progressbar(total as u64);

    for (i,(c, record_dict)) in multi_iter.enumerate(){
        let mut entry: Vec<usize> = Vec::new();
        for s in samplenames.iter(){
            let numi = match record_dict.get(s){
                Some(records) => records.iter()
                    .map(|r|r.UMI)
                    .unique()
                    .count(),
                None => 0
            };
            entry.push(numi)
        }
        result.insert(CB(c), entry);

        if i % 10_000== 0{  // cells iterations are usually rather small, i.e. millions, update more reg
            bar.inc(10_000);
        }
    };
    
    // write to file
    // TODO could be inlined into the above code to instantly write
    let mut fh = File::create(outfile).unwrap();
    let mut header = samplenames.join(",");
    header.push_str(",CB");
    writeln!(fh, "{}", header).unwrap();

    for (cid, numis) in result.iter(){
        // concat with commas
        let mut s = numis.iter().map(|i| i.to_string()).collect::<Vec<String>>().join(",");
        s.push_str(&format!(",{}", cid.0));
        writeln!(fh, "{}", s).unwrap();
    }  
}


pub fn detect_overlap(busfolders: HashMap<String, String>) -> HashMap<Vec<String>, usize> {
    // deprecated
    // mesures the number of ovelapping CUGs across the experiments

    let mut total = 0;
    // TODO: this doesnt check if the EC overlaps
    for v in busfolders.values(){
        println!("determine size of iterator");
        let total_records = BusIteratorBuffered::new(v).groupby_cbumi().count();
        if total< total_records{
            total=total_records
        }
    }
    println!("total records {}", total);

    let multi_iter = CellUmiIteratorMulti::new(&busfolders);
    let bar = get_progressbar(total as u64);
    let mut counter: HashMap<Vec<String>, usize> = HashMap::new();

    for (i,((_cb, _umi), record_dict)) in multi_iter.enumerate(){
        let mut the_set: Vec<String> = record_dict.keys().cloned().collect();
        the_set.sort();
        let val = counter.entry(the_set).or_insert(0);
        *val += 1; 

        if i % 1000000== 0{
            bar.inc(1000000);
        }
    }
    counter
}


/// given a collection of busfiles, FingerPrintHistogram keeps track of the 
/// number of times we see a molecule (CB/UMI/gene) in a particular distribution 
/// across the busfiles. Instead recording individual CUGs, we aggregate all the ones
/// having the same fingerprint, i.e. we create a histogram of fingerprints
/// 
#[derive(Debug)]
pub struct FingerprintHistogram{
    order: Vec<String>,
    histogram: HashMap<Vec<u32>, usize>,
}

#[derive(Eq, PartialEq, Hash, Ord, PartialOrd, Clone, Debug, Copy)]
pub struct AmpFactor(u32);

impl FingerprintHistogram{

    pub fn from_csv(filename: &str) -> Self{
        let fh = File::open(filename).unwrap();
        let mut samplenames: Vec<String> = Vec::new();
        let mut histogram: HashMap<Vec<u32>, usize> = HashMap::new();
        for (i, line) in BufReader::new(fh).lines().flatten().enumerate(){
            if i==0{
                for s in line.split(','){
                    samplenames.push(s.to_string());
                }
                // last one is frequency, remove that
                let fstring = samplenames.pop().unwrap();
                assert_eq!(fstring, "frequency");
            }
            else{
                let mut fingerprint: Vec<u32> = line.split(',').map(|e| e.parse::<u32>().unwrap()).collect();
                // last element is the frequency of the fingerprint
                let freq = fingerprint.pop().unwrap();
                histogram.insert(fingerprint, freq as usize);
            }
        }
        FingerprintHistogram{ 
            order: samplenames,
            histogram,
        }
    }

    pub fn new(sample_order: &[String]) -> Self{
        let hist = HashMap::new();
        FingerprintHistogram{ 
            order: sample_order.to_owned(), // not sure why, linter is suggesting it
            histogram : hist,
        }
    }

    pub fn add(&mut self, record_dict: HashMap<String, Vec<BusRecord>>, ecmapper_dict:&HashMap<String, &Ec2GeneMapper>){
        // creates a fingerprint of the record_dict, updates the counts in the histogram

        let grouped_record_dicts = groupby_gene_across_samples(&record_dict, ecmapper_dict);
        let mut fingerprints: Vec<Vec<u32>> = Vec::new();
        for rd in grouped_record_dicts{

            let fp_hash = make_fingerprint_simple(&rd);
            // turn the hashmap into a vector, sorted acc to order
            let fp:  Vec<_>  = self.order.iter()
                .map(|s| fp_hash.get(s).unwrap_or(&0)).cloned().collect();
            fingerprints.push(fp);
        }

        for fp in fingerprints{
            // update frequency
            let v = self.histogram.entry(fp).or_insert(0);
            *v += 1;            
        }
    }

    pub fn to_csv(&self, outfile: &str){
        let mut fh = File::create(outfile).expect("Cant create file: {outfile}");
        let mut header = self.order.join(",");
        header.push_str(",frequency");
        writeln!(fh, "{}", header).unwrap();

        for (fingerprint, freq) in self.histogram.iter(){
            // concat with commas
            let mut s = fingerprint.iter().map(|i| i.to_string()).collect::<Vec<String>>().join(",");
            s.push_str(&format!(",{}", freq));
            writeln!(fh, "{}", s).unwrap();
        }
    }

    pub fn estimate_sihr(&self) -> f64{
        // number of non-chimeric molecules of amp r
        let mut z_r: HashMap<usize, usize> = HashMap::new();
        
        for (fingerprint, freq ) in self.histogram.iter(){
            let r = fingerprint.iter().map(|x| *x as usize).sum();
            let n_experiments = fingerprint.iter().filter(|x|**x>0).count();
            if n_experiments == 1{
                let v = z_r.entry(r).or_insert(0);
                *v+= freq;
            }
        } 

        let mut m_r: HashMap<usize, usize> = HashMap::new();
        for (fingerprint, freq ) in self.histogram.iter(){
            let r = fingerprint.iter().map(|x| *x as usize).sum();
            let v = m_r.entry(r).or_insert(0);
            *v+= freq;
        } 

        let mut r:Vec<usize> = m_r.iter().map(|(k, _v)| k.to_owned() as usize).collect();
        r.sort();

        let z: Vec<usize> = r.iter().map(|x|*z_r.get(x).unwrap_or(&0)).collect();  // in case there's not a single non-chimeric molecule at amplification r, return 0
        let m: Vec<usize> = r.iter().map(|x|*m_r.get(x).unwrap()).collect();

        let (pmax, _prange, _loglike_range) = phantom_binomial_regression(&z,&m,&r, self.order.len());
        
        let mut fh = File::create("/tmp/phat.csv").unwrap();
        writeln!(fh, "p,logp").unwrap();
        for (p, logp) in izip!(_prange, _loglike_range){
            writeln!(fh, "{},{}", p, logp).unwrap();
        }
        pmax

    }

}


pub fn groupby_gene_across_samples(
    record_dict: &HashMap<String,Vec<BusRecord>>, 
    ecmapper_dict: &HashMap<String, &Ec2GeneMapper>) -> Vec<HashMap<String, BusRecord>> {
    // - build the disjoint set using samplename-> genes, i.e. the disjoint sets elements are samplenames
    // - iterate over the disjoint_set elements (samplenames) and build the Busrecords

    // check if the genes are consistent across samples
    // if so yield the record as is
    // otherwise split into multiple records
    if record_dict.len() == 1{
        // return vec![record_dict];
    };
    
    // give each element in record_dict a unique name ( the key)
    // and store i) the element itself, ii) its genes iii) its samplename
    let mut big_hash:HashMap<String, (&BusRecord, String, &HashSet<u32>)> = HashMap::with_capacity(record_dict.len());
    let mut id_counter: i32 = 0;
    for (sname, records) in record_dict.iter(){
        let ecmapper = ecmapper_dict.get(sname).unwrap();
        for r in records{
            let g = ecmapper.get_genes(r.EC);
            let id_str = id_counter.to_string();
            big_hash.insert(id_str, (r, sname.clone(), g ));
            id_counter+= 1;
        }
    }
    
    // build disjoint set based on samplename and genes
    let mut disjoint_set = DisjointSubsets::new();
    for (id, (_r, _sname, gset)) in big_hash.iter(){
        disjoint_set.add(id.clone(), (*gset).clone());
    }

    // build the emitted dict
    let mut emit_vector: Vec<HashMap<String, BusRecord>> = Vec::new();
    for (ids_of_set_elements, geneset) in disjoint_set.get_disjoint_set_ids_and_set(){
        // these are now all records across the samples that map to the same gene
        // evne if there's multiple records per sample, they can be aggregated into
        // a single record

        // group the records by sample, aggregate
        let the_gene = *geneset.iter().next().unwrap(); // rbitrarily choose the first consissnt gene to mark in busrecords

        let mut sample_grouped: HashMap<String, Vec<BusRecord>> = HashMap::new();
        for el_id in ids_of_set_elements{
            // pop out the element. not needed, but shrinks the map
            let (record, samplename, _genes) = big_hash.remove(&el_id).unwrap();
            
            let rlist = sample_grouped.entry(samplename).or_insert(Vec::new());
            rlist.push(record.clone());
        }

        let mut sample_grouped_aggr: HashMap<String, BusRecord> = HashMap::new();
        for (s, recordlist) in sample_grouped.iter(){
            let freq = recordlist.iter().map(|r| r.COUNT).sum();
            let mut rnew = recordlist.get(0).unwrap().clone();
            rnew.COUNT = freq;
            rnew.EC = the_gene; 
            sample_grouped_aggr.insert(s.to_string(), rnew);
        }
        emit_vector.push(sample_grouped_aggr);
    }
    emit_vector
}


fn create_fingerprint(order: &[String], record_dict: &HashMap<String, Vec<BusRecord>>, ecmapper_dict:&HashMap<String, &Ec2GeneMapper>) -> Vec<Vec<u32>>{

    // turns a record_dict into a fingerprint (an integer vector)
    // a record_Dict can result in mutiple fingerprints if teh genes are inconsistent

    let mut emission: Vec<Vec<u32>> = Vec::new();

    // filter out the cases where a CB/UMI has more than one EC
    let filtered_dict: HashMap<String, BusRecord> = (record_dict).iter()
        .filter(|(_k,v)| v.len()==1)
        .map(|(k, v)| (k.clone(), v[0].clone())) // clone -> pop since we dont use vec after; not critical though
        .collect();

    // for rdict in groupby_gene_simple(filtered_dict, ecmapper_dict){
    for rdict in groupby_gene_even_simpler(filtered_dict, ecmapper_dict){
        let fp_hash = make_fingerprint_simple(&rdict);
        
        // turn the hashmap into a vector, sorted acc to order
        let fp:  Vec<_>  = order.iter()
            .map(|s| fp_hash.get(s).unwrap_or(&0)).cloned().collect();

        emission.push(fp);
    }
    emission
}


fn create_fingerprint2(order: &[String], record_dict: &HashMap<String, Vec<BusRecord>>, ecmapper_dict:&HashMap<String, &Ec2GeneMapper>) -> Vec<Vec<u32>>{

    // turns a record_dict into a fingerprint (an integer vector)
    // a record_Dict can result in mutiple fingerprints if teh genes are inconsistent

    let mut emission: Vec<Vec<u32>> = Vec::new();

    // for rdict in groupby_gene_even_simpler(filtered_dict, ecmapper_dict){
    for rdict in groupby_gene_across_samples(record_dict, ecmapper_dict){

        let fp_hash = make_fingerprint_simple(&rdict);
        
        // turn the hashmap into a vector, sorted acc to order
        let fp:  Vec<_>  = order.iter()
            .map(|s| fp_hash.get(s).unwrap_or(&0)).cloned().collect();

        emission.push(fp);
    }
    emission
}

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

fn map_to_file(h: &HashMap<(AmpFactor, String), f64>, outfile: &str){
    let mut amp: Vec<AmpFactor> = h.keys().map(|(r, _s)| *r).unique().collect();
    let mut samplenames: Vec<String> = h.keys().map(|(_r, s)| s.to_owned()).unique().collect();

    amp.sort();
    samplenames.sort();

    let mut fh = File::create(outfile).unwrap();
    let mut header = "samplename,".to_string();
    let rstring= amp.iter().map(|x| x.0.to_string()).collect::<Vec<String>>().join(",");
    header.push_str(&rstring);
    header.push_str(",frequency");

    writeln!(fh, "{}", header).unwrap();


    for s in samplenames{
        let pi_vec: Vec<String> = amp.iter()
            .map(|r| 
                h.get(&(*r, s.clone())).unwrap().to_string()
            )
            .collect();
        writeln!(fh, "{},{}", s, pi_vec.join(",")).unwrap();
    }
}


#[test]
fn test_posterior(){
    let fph = FingerprintHistogram::from_csv("/home/michi/Dropbox/rustphantompurger/IR56_57_phantom.csv");
    let posterior = PhantomPosterior::new(&fph);
    map_to_file(&posterior.pi_norm, "/tmp/pi_norm.csv");
    map_to_file(&posterior.vr_norm, "/tmp/vr_norm.csv");

    // println!("{:?}", posterior.pi_norm);
}


pub fn make_fingerprint_histogram(busfolders: HashMap<String, BusFolder>) -> FingerprintHistogram{
    // main function here: takes a dict of busfolders, creates fingerprints for each molecule 
    // returns the fingerprints histogram (how often each fingerprint was observed)
    // and the ordering of the fingerprint (i.e. which element corresponds to which experiment)

    // create the EC2gene mappers
    let ecmapper_dict = busfolders.iter()
        .map(|(samplename, bfolder)|
            (samplename.clone(), &bfolder.ec2gene)   //#todo remove clone
        ).collect();

    // a list of busfile-names for the iterator
    let busnames =busfolders.iter()
        .map(|(s, bfolder)| (s.clone(), bfolder.get_busfile()))
        .collect();

    let now = Instant::now();
    let result = _make_fingerprint_histogram(&busnames, &ecmapper_dict);
    let elapsed_time = now.elapsed();
    println!("Ran _make_fingerprint_histogram, took {} seconds.", elapsed_time.as_secs());

    result
}

pub fn _make_fingerprint_histogram(busnames: &HashMap<String, String>, ecmapper_dict: &HashMap<String, &Ec2GeneMapper>) -> FingerprintHistogram{

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

    // the actual workhorse, make_fingerprint_histogram is just a convenient wrapper
    let multi_iter = CellUmiIteratorMulti::new(busnames);
    
    // let bar = ProgressBar::new_spinner();
    // bar.set_style(ProgressStyle::default_bar()
    //     .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos} {per_sec}")
    //     .progress_chars("##-"));
    let bar = get_progressbar(total as u64);

    let mut order: Vec<_> = busnames.keys().cloned().collect();
    order.sort();

    let mut fp_histo = FingerprintHistogram::new(&order);

    for (i,((_cb, _umi), record_dict)) in multi_iter.enumerate(){
        fp_histo.add(record_dict, ecmapper_dict);

        if i % 1000000== 0{
            bar.inc(1000000);
        }
    }
    fp_histo
}

type Fingerprint = HashMap<String, u32>;
fn make_fingerprint(record_dict: &HashMap<String,Vec<BusRecord>>) -> Fingerprint{
    // for a molecule detected over several experiments, 
    // get is fingerprint, i.e. freq across the experiments
    // if there's multiple ECs supporting the read, jsut sum their counts
    let fingerprint: Fingerprint = record_dict
        .iter()
        .map(|(k, v)| {
            (k.clone(), v.iter().map(|r|r.COUNT).sum::<u32>())
        }).collect();
    fingerprint
}

#[inline(never)]
fn make_fingerprint_simple(record_dict: &HashMap<String,BusRecord>) -> Fingerprint{
    // for a molecule detected over several experiments, 
    // get is fingerprint, i.e. freq across the experiments
    let fingerprint: Fingerprint = record_dict
        .iter()
        .map(|(k, v)| {
            (k.clone(), v.COUNT)
        }).collect();
    fingerprint
}

// #[inline(never)]
pub fn groupby_gene_even_simpler(record_dict: HashMap<String,BusRecord>, ecmapper_dict: &HashMap<String, &Ec2GeneMapper>) -> Vec<HashMap<String, BusRecord>> {
    // assuming no multimapped reads, hence record dict contains only single entries (not a list)
    // this uses a much simpler idea without the BusTmp:
    // - build the disjoint set using samplename-> genes, i.e. the disjoint sets elements are samplenames
    // - iterate over the disjoint_set elements (samplenames) and build the Busrecords

    // check if the genes are consistent across samples
    // if so yield the record as is
    // otherwise split into multiple records
    if record_dict.len() == 1{
        return vec![record_dict];
    };
    
    // give each element in record_dict a unique name ( the key)
    // and store i) the element itself, ii) its genes iii) its samplename
    let mut big_hash:HashMap<String, (BusRecord, String, &HashSet<u32>)> = HashMap::with_capacity(record_dict.len());
    for (i, (sname, r)) in record_dict.into_iter().enumerate(){
        let ecmapper = ecmapper_dict.get(&sname).unwrap();
        let g = ecmapper.get_genes(r.EC);
        big_hash.insert(i.to_string(), (r, sname,g ));
    }

    // build disjoint set based on samplename and genes
    let mut disjoint_set = DisjointSubsets::new();
    for (id, (_r, _sname, gset)) in big_hash.iter(){
        disjoint_set.add(id.clone(), (*gset).clone());
    }

    // build the emitted dict
    let mut emit_vector: Vec<HashMap<String, BusRecord>> = Vec::new();
    for ids_of_set_elements in disjoint_set.get_disjoint_set_ids(){
        // these are now all records across the samples that map to the same gene
        // let bus_tuple: Vec<BusTmp> = bus_string_concat.split(SEPARATOR).map(|bstring|BusTmp::parse(bstring).unwrap()).collect();

        let mut emited_dict: HashMap<String, BusRecord> = HashMap::new();
        for el_id in ids_of_set_elements{
            // pop out the element. not needed, but shrinks the map
            let (record, samplename, _genes) = big_hash.remove(&el_id).unwrap();

            if emited_dict.contains_key(&samplename){
                panic!("cant happen, each sample only has one record")
            } 
            else{
                emited_dict.insert(samplename, record);
            }
        }
        emit_vector.push(emited_dict);
    }
    emit_vector
}


#[cfg(test)]
pub mod tests{
    use std::collections::{HashSet, HashMap};
    use rustbustools::{consistent_genes::Ec2GeneMapper, io::{BusRecord, BusFolder, setup_busfile}};
    use statrs::assert_almost_eq;
    use crate::phantompurger::PhantomPosterior;

    use super::{_make_fingerprint_histogram, make_fingerprint_simple, groupby_gene_even_simpler, FingerprintHistogram, groupby_gene_across_samples};
    use super::{make_fingerprint_histogram, detect_cell_overlap};

    #[test]
    fn test_groupby_gene_across_samples(){
        let es1 = create_dummy_ec();
        let es2 = create_dummy_ec();
        let es_dict: HashMap<String, &Ec2GeneMapper> = vec![
            ("s1".to_string(), &es1),
            ("s2".to_string(), &es2),
        ].into_iter().collect();

        // first sample: two records, with consistent gene A
        let r1 =BusRecord{CB: 0, UMI: 1, EC: 0, COUNT: 2, FLAG: 0};
        let r2 =BusRecord{CB: 0, UMI: 1, EC: 2, COUNT: 2, FLAG: 0};

        // second sample: two records, with consistent gene A,B
        let s1 = BusRecord{CB: 0, UMI: 1, EC: 2, COUNT: 3, FLAG: 0};
        let s2 = BusRecord{CB: 0, UMI: 1, EC: 2, COUNT: 3, FLAG: 0};

        let record_dict = vec![
            ("s1".to_string(), vec![r1, r2]),
            ("s2".to_string(), vec![s1, s2]),
        ].into_iter().collect();

        let res = groupby_gene_across_samples(&record_dict, &es_dict);

        // println!("{:?}", res);
        assert_eq!(res.len(), 1);
        // make sure it aggregtes correctly
        let r = res[0].clone();
        assert_eq!(r.get("s1").unwrap().COUNT, 4);
        assert_eq!(r.get("s2").unwrap().COUNT, 6);

        /*
        more complicated example where the ECs dont fully agree
        in sample things indicate A, in sampel B i'ts inconsistent
        looking at all samples together A,B could be an option
         */
        // first sample: two records, with consistent gene A
        let r1 =BusRecord{CB: 0, UMI: 1, EC: 0, COUNT: 2, FLAG: 0}; //A
        let r2 =BusRecord{CB: 0, UMI: 1, EC: 2, COUNT: 2, FLAG: 0}; //A,B

        // second sample: two records, with consistent gene A the other consistent with gene B
        let s1 = BusRecord{CB: 0, UMI: 1, EC: 0, COUNT: 3, FLAG: 0}; // A
        let s2 = BusRecord{CB: 0, UMI: 1, EC: 1, COUNT: 3, FLAG: 0}; //B

        let record_dict = vec![
            ("s1".to_string(), vec![r1, r2]),
            ("s2".to_string(), vec![s1, s2]),
        ].into_iter().collect();
        let res = groupby_gene_across_samples(&record_dict, &es_dict);
        println!("{:?}", res);


    }

    fn create_dummy_ec() ->Ec2GeneMapper{
        let ec0: HashSet<String> = vec!["A".to_string()].into_iter().collect();
        let ec1: HashSet<String> = vec!["B".to_string()].into_iter().collect();
        let ec2: HashSet<String> = vec!["A".to_string(), "B".to_string()].into_iter().collect();
        let ec3: HashSet<String> = vec!["C".to_string(), "D".to_string()].into_iter().collect();

        let ec_dict: HashMap<u32, HashSet<String>> = HashMap::from([
            (0, ec0), // A
            (1, ec1), // B
            (2, ec2), // A,B
            (3, ec3), // C,D
            ]);
        Ec2GeneMapper::new(ec_dict)
    }

    #[test]
    fn test_groupby_gene_simple(){

        // --------------------------------------------
        // two records that share the same EC, should be grouped into a single emission
        let r1 =BusRecord{CB: 0, UMI: 1, EC: 0, COUNT: 2, FLAG: 0};
        let s1 = BusRecord{CB: 0, UMI: 1, EC: 0, COUNT: 3, FLAG: 0};

        let es1 = create_dummy_ec();
        let es2 = create_dummy_ec();
        let es3 = create_dummy_ec();
        let es_dict: HashMap<String, &Ec2GeneMapper> = vec![
            ("s1".to_string(), &es1),
            ("s2".to_string(), &es2),
            ("s3".to_string(), &es3)
            ]
        .into_iter().collect();


        let record_dict = vec![
            ("s1".to_string(), r1),
            ("s2".to_string(), s1),
        ].into_iter().collect();
        // let res = groupby_gene_simple(record_dict, &es_dict);
        let res = groupby_gene_even_simpler(record_dict, &es_dict);
        println!("{:?}", res);
        assert_eq!(res.len(), 1);
        assert_eq!(res[0].len(), 2);

        // --------------------------------------------
        // two records, different EC, but same gene, should be grouped into a single emission
        let r1 =BusRecord{CB: 0, UMI: 1, EC: 0, COUNT: 1, FLAG: 0};
        let s1 = BusRecord{CB: 0, UMI: 1, EC: 2, COUNT: 1, FLAG: 0};
        let record_dict = vec![
            ("s1".to_string(), r1),
            ("s2".to_string(), s1),
        ].into_iter().collect();
        // let res = groupby_gene_simple(record_dict, &es_dict);
        let res = groupby_gene_even_simpler(record_dict, &es_dict);

        println!("{:?}", res);
        assert_eq!(res.len(), 1);
        assert_eq!(res[0].len(), 2);

        // --------------------------------------------
        // two records, different EC, and inconsistnet gene, should be grouped into a TWO emissions
        let r1 =BusRecord{CB: 0, UMI: 1, EC: 0, COUNT: 1, FLAG: 0};
        let s1 = BusRecord{CB: 0, UMI: 1, EC: 1, COUNT: 1, FLAG: 0};
        let record_dict = vec![
            ("s1".to_string(), r1),
            ("s2".to_string(), s1),
        ].into_iter().collect();
        // let res = groupby_gene_simple(record_dict, &es_dict);
        let res = groupby_gene_even_simpler(record_dict, &es_dict);

        println!("{:?}", res);
        assert_eq!(res.len(), 2);
        assert_eq!(res[0].len(), 1);
        assert_eq!(res[1].len(), 1);

        // --------------------------------------------
        // three records, A, B, (A,B). should be yielded as a single emssion
        let r1 =BusRecord{CB: 0, UMI: 1, EC: 0, COUNT: 1, FLAG: 0};
        let s1 = BusRecord{CB: 0, UMI: 1, EC: 1, COUNT: 1, FLAG: 0};
        let t1 = BusRecord{CB: 0, UMI: 1, EC: 2, COUNT: 1, FLAG: 0};
        let record_dict = vec![
            ("s1".to_string(), r1),
            ("s2".to_string(), s1),
            ("s3".to_string(), t1),
        ].into_iter().collect();
        // let res = groupby_gene_simple(record_dict, &es_dict);
        let res = groupby_gene_even_simpler(record_dict, &es_dict);

        println!("{:?}", res);
        assert_eq!(res.len(), 1);
        assert_eq!(res[0].len(), 3);

    }

    #[test]
    fn test_make_fingerprint_histogram(){
        use rustbustools::io::setup_busfile;

        // a pair, same EC
        let r1 =BusRecord{CB: 0, UMI: 1, EC: 0, COUNT: 2, FLAG: 0};
        let s1 = BusRecord{CB: 0, UMI: 1, EC: 0, COUNT: 3, FLAG: 0};
        
        // singleton in sample1
        let r2 =BusRecord{CB: 1, UMI: 2, EC: 0, COUNT: 12, FLAG: 0}; 
        
        // singleton in sample2
        let s2 = BusRecord{CB: 1, UMI: 2, EC: 1, COUNT: 10, FLAG: 0};


        // pair but only gene overlap
        let r3 =BusRecord{CB: 1, UMI: 3, EC: 0, COUNT:  2, FLAG: 0}; 
        let s3 = BusRecord{CB: 1, UMI: 3, EC: 2, COUNT:  3, FLAG: 0}; 

        //overall we should get 
        // [12, 0] = 1
        // [0, 10] = 1
        // [2, 3] = 2

        let v1 = vec![r1.clone(),r2.clone(),r3.clone()];
        let v2 = vec![s1.clone(),s2.clone(),s3.clone()];

        // write the records to file
        let (busname1, _dir1) =setup_busfile(&v1);
        let (busname2, _dir2) =setup_busfile(&v2);

        let hashmap = HashMap::from([
            ("s1".to_string(), busname1.to_string()),
            ("s2".to_string(), busname2.to_string())
        ]);

        let es1 = create_dummy_ec();
        let es2 = create_dummy_ec();
        let es_dict: HashMap<String, &Ec2GeneMapper> = vec![
            ("s1".to_string(), &es1),
            ("s2".to_string(), &es2)]
            .into_iter().collect();

        let s = _make_fingerprint_histogram(&hashmap, &es_dict);
        println!("{:?}", s);

        let e = *s.histogram.get(&vec![12, 0]).unwrap();
        assert_eq!(e, 1);

        let e = *s.histogram.get(&vec![0, 10]).unwrap();
        assert_eq!(e, 1);

        let e = *s.histogram.get(&vec![2, 3]).unwrap();
        assert_eq!(e, 2);

        s.to_csv("/tmp/finger.csv")

        // let s = detect_overlap(hashmap);
    }

    #[test]
    fn test_make_fingerprint_simple(){

        let r1 =BusRecord{CB: 0, UMI: 1, EC: 0, COUNT: 2, FLAG: 0};
        let s1 = BusRecord{CB: 0, UMI: 1, EC: 0, COUNT: 3, FLAG: 0};

        let record_dict = vec![
            ("s1".to_string(), r1),
            ("s2".to_string(), s1),
        ].into_iter().collect();
        let res = make_fingerprint_simple(&record_dict);
        println!("{:?}", res);

        let exp: HashMap<_,_> = vec![("s1".to_string(), 2), ("s2".to_string(), 3)].into_iter().collect();
        assert_eq!(res, exp);
    }


    // #[test]
    pub fn testing2(){
        
        let t2g = "/home/michi/bus_testing/transcripts_to_genes.txt";
        // let b1 = BusFolder::new("/home/michi/bus_testing/bus_output/", t2g);
        let b1 = BusFolder::new("/home/michi/bus_testing/bus_output_short/", t2g);
        let b2 = BusFolder::new("/home/michi/bus_testing/bus_output_short/", t2g);

        // let hashmap = HashMap::from([
        //     ("full".to_string(), "/home/michi/bus_testing/bus_output/output.corrected.sort.bus".to_string()),
        //     ("short".to_string(), "/home/michi/bus_testing/bus_output_short/output.corrected.sort.bus".to_string())
        // ]);
        let hashmap = HashMap::from([
            ("full".to_string(), b1),
            ("short".to_string(), b2)
        ]);

        let s = make_fingerprint_histogram(hashmap);

        s.to_csv("/tmp/testing2.csv");
        println!("{:?}", s)
    }

    // #[test]
    pub fn test_detect_cell_overlap(){
        // let t2g = "/home/michi/bus_testing/transcripts_to_genes.txt";
        // let b1 = BusFolder::new("/home/michi/bus_testing/bus_output/", t2g);
        // let b1 = BusFolder::new("/home/michi/bus_testing/bus_output_short/", t2g);
        // let b2 = BusFolder::new("/home/michi/bus_testing/bus_output_short/", t2g);
        // let busfolders = HashMap::from([
        //     ("full".to_string(), b1),
        //     ("short".to_string(), b2)
        // ]);

        let busfolders = HashMap::from([
            ("full".to_string(), "/home/michi/bus_testing/bus_output/output.corrected.sort.bus".to_string()),
            ("short".to_string(), "/home/michi/bus_testing/bus_output_short/output.corrected.sort.bus".to_string())
        ]);

        detect_cell_overlap(busfolders, "/tmp/test_detect_cell_overlap.csv")
    }

    #[test]
    pub fn test_csv_read_write(){

        let order = vec!["s1", "s2"].into_iter().map(|x| x.to_string()).collect();
        let mut histogram: HashMap<Vec<u32>, usize> = HashMap::new();

        histogram.insert(
            vec![1, 0],
            1
        );
        histogram.insert(
            vec![0, 1],
            1
        );
        histogram.insert(
            vec![1, 1],
            1
        );

        let fph = FingerprintHistogram {
            order,
            histogram,
        };

        fph.to_csv("/tmp/finger.csv");

        println!("{:?}", fph);
        let t = FingerprintHistogram::from_csv("/tmp/finger.csv");
        println!("{:?}", t);
        
        assert_eq!(fph.histogram, t.histogram);

        let p = t.estimate_sihr();
        println!("{}", p);
    }

    #[test]
    pub fn test_estimate_sihr(){

        let order = vec!["s1", "s2"].into_iter().map(|x| x.to_string()).collect();
        let mut histogram: HashMap<Vec<u32>, usize> = HashMap::new();

        histogram.insert(
            vec![1, 0],
            1
        );
        histogram.insert(
            vec![0, 1],
            1
        );
        histogram.insert(
            vec![1, 1],
            1
        );

        let fph = FingerprintHistogram {
            order,
            histogram,
        };

        let p = fph.estimate_sihr();
        println!("{}", p);

        assert_almost_eq!(p,  0.5, 0.001);

    }

    #[test]
    fn test_filter(){

        // exclusive in 1
        let r1 = BusRecord{CB: 0, UMI: 0, EC: 0, COUNT: 100, FLAG: 0};
        // exclusuve in 2
        let s1 = BusRecord{CB: 0, UMI: 21, EC: 1, COUNT: 100, FLAG: 0};
        // exclusuve in 2
        let t1 = BusRecord{CB: 0, UMI: 1, EC: 0, COUNT: 100, FLAG: 0};

        // shared [10,1] in 1/2
        let r2 = BusRecord{CB: 1, UMI: 0, EC: 0, COUNT: 100, FLAG: 0};
        let s2 = BusRecord{CB: 1, UMI: 0, EC: 2, COUNT: 1, FLAG: 0};
        let t2 = BusRecord{CB: 1, UMI: 0, EC: 2, COUNT: 1, FLAG: 0};

        let (b1, tdir1) = setup_busfile(&vec![r1,r2]);
        let (b2, tdir2) = setup_busfile(&vec![s1,s2]);
        let (b3, tdir3) = setup_busfile(&vec![t1,t2]);

        println!("{}", b1);

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

}
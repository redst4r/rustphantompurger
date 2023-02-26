// use core::panic;
use std::{collections::{HashMap, HashSet}, fs::File, io::{Write, BufReader, BufRead}, time::Instant, hash::Hash};
use crate::{disjoint::{DisjointSubsets}, binomialreg::phantom_binomial_regression};
use indicatif::{ProgressBar, ProgressStyle};
use itertools::izip;
use rustbustools::{bus_multi::{CellUmiIteratorMulti, CellIteratorMulti}, io::{BusRecord, BusFolder, BusWriter}, iterators::CbUmiIterator, utils::{get_progressbar, argsort::argmax_float}, consistent_genes::Ec2GeneMapper};


#[derive(Eq, PartialEq, Hash, Ord, PartialOrd, Clone, Debug,Copy)]
pub struct CB(u64);

pub fn detect_cell_overlap(busfolders: HashMap<String, String>, outfile: &str) {
    // for each CB calcualte the number of UMIs per experiment/busfile
    // if extensive swapping has occured, we'll see shared CBs across experiments
    // with correlated #umis

    let samplenames: Vec<String> = busfolders.keys().cloned().collect();
    let multi_iter = CellIteratorMulti::new(&busfolders);
    let mut result: HashMap<CB, Vec<usize>> = HashMap::new();

    for (c, record_dict) in multi_iter{
        let mut entry: Vec<usize> = Vec::new();
        for s in samplenames.iter(){
            let numi = match record_dict.get(s){
                Some(records) => records.len(),
                None => 0
            };
            entry.push(numi)
        }
        result.insert(CB(c), entry);
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
        let cbumi_iter_tmp = CbUmiIterator::new(v);
        println!("determine size of iterator");
        // let now = Instant::now();
        let total_records = cbumi_iter_tmp.count();
        if total< total_records{
            total=total_records
        }
    }
    println!("total records {}", total);


    let multi_iter = CellUmiIteratorMulti::new(&busfolders);
    
    let bar = get_progressbar(total as u64);
    
    // let mut counter: HashMap<HashSet<String>, usize> = HashMap::new();
    let mut counter: HashMap<Vec<String>, usize> = HashMap::new();


    for (i,((_cb, _umi), record_dict)) in multi_iter.enumerate(){
        // let the_Set: HashSet<String> = record_dict.keys().cloned().collect();
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
    drop_multimapped: bool
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
            drop_multimapped: true
        }

    }

    pub fn new(sample_order: &[String]) -> Self{
        let hist = HashMap::new();
        FingerprintHistogram{ 
            order: sample_order.to_owned(), // not sure why, linter is suggesting it
            histogram : hist,
            drop_multimapped: true
        }
    }

    pub fn add(&mut self, record_dict: HashMap<String, Vec<BusRecord>>, ecmapper_dict:&HashMap<String, &Ec2GeneMapper>){
        // creates a fingerprint of the record_dict, updates the counts in the histogram

        let fingerprints = create_fingerprint(&self.order, &record_dict, ecmapper_dict);

        for fp in fingerprints{
            // update frequency
            let v = self.histogram.entry(fp).or_insert(0);
            *v += 1;            
        }
    }

    pub fn to_csv(&self, outfile: &str){
        let mut fh = File::create(outfile).unwrap();
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

        println!("r: {:?}", r);
        println!("zr: {:?}", z_r);
        println!("mr: {:?}", m_r);


        let z: Vec<usize> = r.iter().map(|x|*z_r.get(x).unwrap_or(&0)).collect();  // in case there's not a single non-chimeric molecule at amplification r, return 0
        let m: Vec<usize> = r.iter().map(|x|*m_r.get(x).unwrap()).collect();

        println!("z: {:?}", z);
        println!("m: {:?}", m);

        let (pmax, _prange, _loglike_range) = phantom_binomial_regression(&z,&m,&r, self.order.len());
        1.0 - pmax

    }

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

pub struct PhantomPosterior{
    order: Vec<String>,
    // posterior_cache: HashMap<Vec<u32>, Vec<f64>>,  // for each fingerprint a vector of posterior probs
    pi_norm: HashMap<(AmpFactor, String), f64>,
    p_no_hop: f64
}

impl PhantomPosterior{

    pub fn new(fph: &FingerprintHistogram) -> Self{

        let order: Vec<String> = fph.order.iter().cloned().collect();
        let n_samples = order.len() as f64;
        let p = 1.0-fph.estimate_sihr(); // the prop of not hopping

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

        // convert to pi_rs
        let pi_norm: HashMap<(AmpFactor, String), f64> = vr_norm
            .into_iter()
            .map(|((r, s), f)|
                ((r, s), (f * (n_samples-1.0) + (p-1.0))/ (n_samples*p -1.0))
            )
            .collect();        

        // let posterior_cache = HashMap::new();
        PhantomPosterior{
            order,
            pi_norm,  // for each fingerprint a vector of posterior probs
            // posterior_cache,
            p_no_hop: p
        }
    }

    fn posterior_internal(&self, fingerprint: &Vec<u32>) -> Vec<f64>{
        let n_samples = self.order.len() as f64;

        assert!(fingerprint.len()==self.order.len());

        let r = AmpFactor(fingerprint.iter().sum());

        let mut posterior: Vec<f64> = Vec::with_capacity(self.order.len());

        for (sname, y) in izip![self.order.iter(), fingerprint.iter() ]{

            let pi = *self.pi_norm.get(&(r, sname.to_string())).unwrap();

            let base = (n_samples-1.0)/((1.0/ self.p_no_hop) -1.0);
            let post = base.powi(*y as i32) * pi;
            posterior.push(post);
        }
        let norm_constant: f64 = posterior.iter().sum();
        let posterior_normed: Vec<f64> = posterior.into_iter().map(|v| (v/norm_constant)).collect();
        posterior_normed
    }

    pub fn posterior(&self, fingerprint: Fingerprint) -> HashMap<String, f64>{

        let mut fp_numeric: Vec<u32> = Vec::with_capacity(self.order.len());
        for s in self.order.iter(){
            let y = match fingerprint.get(s){
                Some(x) => x,
                None => &0
            };
            fp_numeric.push(*y);
        }
        let posterior_numeric = self.posterior_internal(&fp_numeric);

        let mut posterior_map: HashMap<String, f64> = HashMap::new();

        for (sname,post) in izip![self.order.iter(), posterior_numeric.iter()]{
            posterior_map.insert(sname.clone(), *post);
        }
        posterior_map
    }


    pub fn filter_busfiles(&self, input_busfolders: HashMap<String, BusFolder>, output_busfolders: HashMap<String, String>){
        // for each fingerprint a vector of posterior probs
        let mut posterior_cache: HashMap<Vec<u32>, Vec<f64>> = HashMap::new();  

        // create the EC2gene mappers
        let ecmapper_dict = input_busfolders.iter()
        .map(|(samplename, bfolder)|
            (samplename.clone(), &bfolder.ec2gene)   //#todo remove clone
        ).collect();

        // a list of busfile-names for the iterator
        let busnames =input_busfolders.iter()
            .map(|(s, bfolder)| (s.clone(), bfolder.get_busfile()))
            .collect();
    
        let mut buswriters: HashMap<String,BusWriter> = output_busfolders.iter()
            .map(|(sname, fname)| 
                (
                    sname.to_owned(), 
                    BusWriter::new(
                        &fname, 
                        input_busfolders.get(sname).unwrap().get_bus_header()
                    )
                ) 
            )
            .collect();

        let multi_iter = CellUmiIteratorMulti::new(&busnames);
        let mut records_resolved = 0;
        let mut records_total = 0;

        for (_i,(_cb_umi, record_dict)) in multi_iter.enumerate(){
            records_total += 1;

            let fps = create_fingerprint(&self.order, &record_dict, &ecmapper_dict);

            // if we get a record dict that turns into more than one fingerprint, just ignore
            if fps.len() >1 {
                continue
            }
            let fp = fps.first().unwrap();


            if !posterior_cache.contains_key(fp){
                let pvec = self.posterior_internal(fp);
                posterior_cache.insert(fp.clone(), pvec);
            }
            let posterior_vec = posterior_cache.get(fp).unwrap();


            let (ix, pmax) = argmax_float(posterior_vec);
            let sample_max = self.order[ix].clone();
            // if we find an unambiguous assignment, write the molecule to that sample
            // otherwise drop the molecule in all samples
            if pmax > 0.999999{
                let wr = buswriters.get_mut(&sample_max).unwrap();
                let r = record_dict.get(&sample_max).unwrap();
                wr.write_records(r);
                records_resolved+=1;
            }
        }
        println!("{records_resolved}/{records_total} records corrected/written")
    }
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

    // the actual workhorse, make_fingerprint_histogram is just a convenient wrapper
    let multi_iter = CellUmiIteratorMulti::new(busnames);
    
    let bar = ProgressBar::new_spinner();
    bar.set_style(ProgressStyle::default_bar()
        .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos} {per_sec}")
        .progress_chars("##-"));
        
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
    use rustbustools::{consistent_genes::Ec2GeneMapper, io::{BusRecord, BusFolder}};
    use statrs::assert_almost_eq;
    use super::{_make_fingerprint_histogram, make_fingerprint_simple, groupby_gene_even_simpler, FingerprintHistogram};
    use super::{make_fingerprint_histogram, detect_cell_overlap};

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
            drop_multimapped: true,
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
            drop_multimapped: true,
        };

        let p = fph.estimate_sihr();
        println!("{}", p);

        assert_almost_eq!(p,  0.5, 0.001);

    }
}
use crate::disjoint::DisjointSubsets;
use crate::utils::ec_mapper_dict_from_busfolders;
use crate::binomialreg::phantom_binomial_regression;
use bustools::consistent_genes::Ec2GeneMapper;
use bustools::io::{BusFolder, BusRecord};
use bustools::merger::MultiIterator;
use bustools::utils::get_spinner;
use bustools::{
    consistent_genes::{GeneId, Genename, EC},
    iterators::CbUmiGroupIterator,
};
use serde::{Deserialize, Serialize};
use std::{
    collections::{HashMap, HashSet},
    fs::File,
    hash::Hash,
    io::{BufRead, BufReader, Write},
    time::Instant,
};

/// given a collection of busfiles, FingerPrintHistogram keeps track of the
/// number of times we see a molecule (CB/UMI/gene) in a particular distribution
/// across the busfiles. Instead recording individual CUGs, we aggregate all the ones
/// having the same fingerprint across samples, i.e. we create a histogram of fingerprints
#[derive(Debug)]
pub struct FingerprintHistogram {
    // Names of the samples across which we calcualte the fingerprints
    pub(crate) samplenames: Vec<String>,
    // Histogram FingerPrint -> Frequency
    pub(crate) histogram: HashMap<Vec<u32>, usize>,
}

#[derive(Eq, PartialEq, Hash, Ord, PartialOrd, Clone, Debug, Copy, Deserialize, Serialize)]
pub struct AmpFactor(u32);

impl FingerprintHistogram {

    pub fn get_samplenames(&self) -> Vec<String> {
        self.samplenames.clone()
    }

    pub fn get_histogram(&self) -> HashMap<Vec<u32>, usize> {
        self.histogram.clone()
    }
        
    /// load histogram from disk
    pub fn from_csv(filename: &str) -> Self {
        let fh = File::open(filename).unwrap();
        let mut samplenames: Vec<String> = Vec::new();
        let mut histogram: HashMap<Vec<u32>, usize> = HashMap::new();
        for (i, line) in BufReader::new(fh).lines().flatten().enumerate() {
            if i == 0 {
                for s in line.split(',') {
                    samplenames.push(s.to_string());
                }
                // last one is frequency, remove that
                let fstring = samplenames.pop().unwrap();
                assert_eq!(fstring, "frequency");
            } else {
                let mut fingerprint: Vec<u32> =
                    line.split(',').map(|e| e.parse::<u32>().unwrap()).collect();
                // last element is the frequency of the fingerprint
                let freq = fingerprint.pop().unwrap();
                histogram.insert(fingerprint, freq as usize);
            }
        }
        FingerprintHistogram {
            samplenames,
            histogram,
        }
    }

    pub fn new(sample_order: &[String]) -> Self {
        let hist = HashMap::new();
        FingerprintHistogram {
            samplenames: sample_order.to_owned(), // not sure why, linter is suggesting it
            histogram: hist,
        }
    }

    /// adds a CUG (observed over several busfiles) into the histogram
    pub fn add(
        &mut self,
        record_dict: HashMap<String, Vec<BusRecord>>,
        ecmapper_dict: &HashMap<String, &Ec2GeneMapper>,
    ) {
        // creates a fingerprint of the record_dict, updates the counts in the histogram

        let grouped_record_dicts = groupby_gene_across_samples(&record_dict, ecmapper_dict);
        let mut fingerprints: Vec<Vec<u32>> = Vec::new();
        for rd in grouped_record_dicts {
            let fp_hash = make_fingerprint_simple(&rd);
            // turn the hashmap into a vector, sorted acc to order
            let fp: Vec<_> = self
                .samplenames
                .iter()
                .map(|s| fp_hash.get(s).unwrap_or(&0))
                .cloned()
                .collect();
            fingerprints.push(fp);
        }

        for fp in fingerprints {
            // update frequency
            let v = self.histogram.entry(fp).or_insert(0);
            *v += 1;
        }
    }

    /// save the histogram to disk
    pub fn to_csv(&self, outfile: &str) {
        let mut fh = File::create(outfile).expect("Cant create file: {outfile}");
        let mut header = self.samplenames.join(",");
        header.push_str(",frequency");
        writeln!(fh, "{}", header).unwrap();

        for (fingerprint, freq) in self.histogram.iter() {
            // concat with commas
            let mut s = fingerprint
                .iter()
                .map(|i| i.to_string())
                .collect::<Vec<String>>()
                .join(",");
            s.push_str(&format!(",{}", freq));
            writeln!(fh, "{}", s).unwrap();
        }
    }


    /// get the number of non-chimeric molecules as a function of r
    fn get_z_r(&self) -> HashMap<AmpFactor, usize>{
        let mut z_r: HashMap<AmpFactor, usize> = HashMap::new();

        for (fingerprint, freq) in self.histogram.iter() {
            let r: usize = fingerprint.iter().map(|x| *x as usize).sum();
            let n_experiments = fingerprint.iter().filter(|x| **x > 0).count();
            if n_experiments == 1 {
                let v = z_r.entry(AmpFactor(r as u32)).or_insert(0);
                *v += freq;
            }
        }      
        z_r  
    }

    /// get the number of total molecules as a function of r
    fn get_m_r(&self) -> HashMap<AmpFactor, usize>{

        let mut m_r: HashMap<AmpFactor, usize> = HashMap::new();
        for (fingerprint, freq) in self.histogram.iter() {
            let r: usize = fingerprint.iter().map(|x| *x as usize).sum();
            let v = m_r.entry(AmpFactor(r as u32)).or_insert(0);
            *v += freq;
        }
        m_r
    }
    
    /// Estimate SIHR from the histogram via binomial regression
    /// TODO: the paper suggests to filter some values:
    /// - remove samples (rows over r) with less than 10 observations (m_r <10)
    /// - remove samples with r>25
    pub fn estimate_sihr(&self) -> f64 {
        // number of non-chimeric molecules of amp r

        let z_r = self.get_z_r();
        let m_r = self.get_m_r();

        let mut r: Vec<AmpFactor> = m_r.keys().map(|k| k.to_owned()).collect();
        r.sort();

        let z: Vec<usize> = r.iter().map(|x| *z_r.get(x).unwrap_or(&0)).collect(); // in case there's not a single non-chimeric molecule at amplification r, return 0
        let m: Vec<usize> = r.iter().map(|x| *m_r.get(x).unwrap()).collect();

        let r_usize: Vec<_> = r.iter().map(|x| x.0 as usize).collect(); // convert from AmpFactor -> usize
        let (pmax, _prange, _loglike_range) =
            phantom_binomial_regression(&z, &m, &r_usize, self.samplenames.len());

        // let mut fh = File::create("/tmp/phat.csv").unwrap();
        // writeln!(fh, "p,logp").unwrap();
        // for (p, logp) in izip!(_prange, _loglike_range) {
        //     writeln!(fh, "{},{}", p, logp).unwrap();
        // }
        pmax
    }
}

pub fn groupby_gene_across_samples(
    record_dict: &HashMap<String, Vec<BusRecord>>,
    ecmapper_dict: &HashMap<String, &Ec2GeneMapper>,
) -> Vec<HashMap<String, BusRecord>> {
    // - build the disjoint set using samplename-> genes, i.e. the disjoint sets elements are samplenames
    // - iterate over the disjoint_set elements (samplenames) and build the Busrecords

    // check if the genes are consistent across samples
    // if so yield the record as is
    // otherwise split into multiple records
    if record_dict.len() == 1 {
        // return vec![record_dict];
    };

    // give each element in record_dict a unique name ( the key)
    // and store i) the element itself, ii) its genes iii) its samplename
    let mut big_hash: HashMap<String, (&BusRecord, String, &HashSet<GeneId>)> =
        HashMap::with_capacity(record_dict.len());
    let mut id_counter: i32 = 0;
    for (sname, records) in record_dict.iter() {
        let ecmapper = ecmapper_dict.get(sname).unwrap();
        for r in records {
            let g = ecmapper.get_genes(EC(r.EC));
            let id_str = id_counter.to_string();
            big_hash.insert(id_str, (r, sname.clone(), g));
            id_counter += 1;
        }
    }

    // build disjoint set based on samplename and genes
    let mut disjoint_set = DisjointSubsets::new();
    for (id, (_r, _sname, gset)) in big_hash.iter() {
        disjoint_set.add(id.clone(), (*gset).clone());
    }

    // build the emitted dict
    let mut emit_vector: Vec<HashMap<String, BusRecord>> = Vec::new();
    for (ids_of_set_elements, _geneset) in disjoint_set.get_disjoint_set_ids_and_set() {
        // these are now all records across the samples that map to the same gene
        // evne if there's multiple records per sample, they can be aggregated into
        // a single record

        // group the records by sample, aggregate

        // deprecated:
        // let the_gene = *geneset.iter().next().unwrap(); // rbitrarily choose the first consissnt gene to mark in busrecords

        let mut sample_grouped: HashMap<String, Vec<BusRecord>> = HashMap::new();
        for el_id in ids_of_set_elements {
            // pop out the element. not needed, but shrinks the map
            let (record, samplename, _genes) = big_hash.remove(&el_id).unwrap();

            let rlist = sample_grouped.entry(samplename).or_default();
            rlist.push(record.clone());
        }

        let mut sample_grouped_aggr: HashMap<String, BusRecord> = HashMap::new();
        for (s, recordlist) in sample_grouped.iter() {
            let freq = recordlist.iter().map(|r| r.COUNT).sum();
            let mut rnew = recordlist.first().unwrap().clone();
            rnew.COUNT = freq;

            // bad idea here: this turns the BusRecord into something not quantifiable by bustools!!
            // rnew.EC = the_gene;
            // rather just leave the EC as is. It's guaranteed to be consistent eith the_gene
            // but might not be unique to the_gene
            // TODO: figure out a EC that maps ONLY to the_gene
            sample_grouped_aggr.insert(s.to_string(), rnew);
        }
        emit_vector.push(sample_grouped_aggr);
    }
    emit_vector
}



/// main function here: takes a dict of busfolders, creates fingerprints for each molecule
/// returns the fingerprints histogram (how often each fingerprint was observed)
/// and the ordering of the fingerprint (i.e. which element corresponds to which experiment)
pub fn make_fingerprint_histogram(busfolders: HashMap<String, BusFolder>, t2g_file: &str) -> FingerprintHistogram {

    // create the EC2gene mappers
    // silly, cant create one where ECMapper is a reference
    // need to instantiate/own it, then create ref

    println!("Making EC mappers");
    let ec_tmp = ec_mapper_dict_from_busfolders(&busfolders, t2g_file );
    let ecmapper_dict: HashMap<String, &Ec2GeneMapper> = ec_tmp.iter().map(|(name, d)| (name.clone(), d )).collect();

    let now = Instant::now();
    let result = _make_fingerprint_histogram(&busfolders, &ecmapper_dict);
    let elapsed_time = now.elapsed();
    println!(
        "Ran _make_fingerprint_histogram, took {} seconds.",
        elapsed_time.as_secs()
    );
    result
}

fn _make_fingerprint_histogram(
    busfolders: &HashMap<String, BusFolder>,
    ecmapper_dict: &HashMap<String, &Ec2GeneMapper>,
) -> FingerprintHistogram {

    // the actual workhorse, make_fingerprint_histogram is just a convenient wrapper
    let iterators = busfolders.iter().map(|(k,v)| (k.clone(), v.get_iterator().groupby_cbumi())).collect();
    let multi_iter = MultiIterator::new(iterators);

    // let bar = get_progressbar(total as u64);
    let bar = get_spinner();
    
    let mut order: Vec<_> = busfolders.keys().cloned().collect();
    order.sort();

    let mut fp_histo = FingerprintHistogram::new(&order);

    for (i, ((_cb, _umi), record_dict)) in multi_iter.enumerate() {
        fp_histo.add(record_dict, ecmapper_dict);

        if i % 1000000 == 0 {
            bar.inc(1000000);
        }
    }
    fp_histo
}

pub type Fingerprint = HashMap<String, u32>;

//#[inline(never)]
pub fn make_fingerprint_simple(record_dict: &HashMap<String, BusRecord>) -> Fingerprint {
    // for a molecule detected over several experiments,
    // get is fingerprint, i.e. freq across the experiments
    let fingerprint: Fingerprint = record_dict
        .iter()
        .map(|(k, v)| (k.clone(), v.COUNT))
        .collect();
    fingerprint
}

pub fn create_dummy_ec() -> Ec2GeneMapper {
    let ec0: HashSet<Genename> = vec![Genename("A".to_string())].into_iter().collect();
    let ec1: HashSet<Genename> = vec![Genename("B".to_string())].into_iter().collect();
    let ec2: HashSet<Genename> = vec![Genename("A".to_string()), Genename("B".to_string())].into_iter().collect();
    let ec3: HashSet<Genename> = vec![Genename("C".to_string()), Genename("D".to_string())].into_iter().collect();

    let ec_dict: HashMap<EC, HashSet<Genename>> = HashMap::from([
        (EC(0), ec0), // A
        (EC(1), ec1), // B
        (EC(2), ec2), // A,B
        (EC(3), ec3), // C,D
    ]);
    Ec2GeneMapper::new(ec_dict)
}

#[cfg(test)]
pub mod tests {
    use crate::phantompurger::create_dummy_ec;
    use bustools::{
        consistent_genes::Ec2GeneMapper,
        io::{BusFolder, BusRecord},
    };
    use statrs::assert_almost_eq;
    use std::collections::HashMap;

    use super::{
        _make_fingerprint_histogram, groupby_gene_across_samples, make_fingerprint_simple,
        FingerprintHistogram,
    };
    use super::make_fingerprint_histogram;

    #[test]
    fn test_groupby_gene_across_samples() {
        let es1 = create_dummy_ec();
        let es2 = create_dummy_ec();
        let es_dict: HashMap<String, &Ec2GeneMapper> =
            vec![("s1".to_string(), &es1), ("s2".to_string(), &es2)]
                .into_iter()
                .collect();

        // first sample: two records, with consistent gene A
        let r1 =BusRecord{CB: 0, UMI: 1, EC: 0, COUNT: 2, FLAG: 0};
        let r2 =BusRecord{CB: 0, UMI: 1, EC: 2, COUNT: 2, FLAG: 0};

        // second sample: two records, with consistent gene A,B
        let s1 = BusRecord{CB: 0, UMI: 1, EC: 2, COUNT: 3, FLAG: 0};
        let s2 = BusRecord{CB: 0, UMI: 1, EC: 2, COUNT: 3, FLAG: 0};

        let record_dict = vec![
            ("s1".to_string(), vec![r1, r2]),
            ("s2".to_string(), vec![s1, s2]),
        ]
        .into_iter()
        .collect();

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
        ]
        .into_iter()
        .collect();
        let res = groupby_gene_across_samples(&record_dict, &es_dict);
        println!("{:?}", res);
    }

    #[test]
    fn test_make_fingerprint_histogram() {
        use bustools::io::setup_busfile;

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

        let v1 = vec![r1.clone(), r2.clone(), r3.clone()];
        let v2 = vec![s1.clone(), s2.clone(), s3.clone()];

        // write the records to file
        let (busname1, _dir1) = setup_busfile(&v1);
        let (busname2, _dir2) = setup_busfile(&v2);

        let hashmap = HashMap::from([
            ("s1".to_string(), BusFolder::new(&busname1)),
            ("s2".to_string(), BusFolder::new(&busname2)),
        ]);

        let es1 = create_dummy_ec();
        let es2 = create_dummy_ec();
        let es_dict: HashMap<String, &Ec2GeneMapper> =
            vec![("s1".to_string(), &es1), ("s2".to_string(), &es2)]
                .into_iter()
                .collect();

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
    fn test_make_fingerprint_simple() {

        let r1 =BusRecord{CB: 0, UMI: 1, EC: 0, COUNT: 2, FLAG: 0};
        let s1 = BusRecord{CB: 0, UMI: 1, EC: 0, COUNT: 3, FLAG: 0};

        let record_dict = vec![("s1".to_string(), r1), ("s2".to_string(), s1)]
            .into_iter()
            .collect();
        let res = make_fingerprint_simple(&record_dict);
        println!("{:?}", res);

        let exp: HashMap<_, _> = vec![("s1".to_string(), 2), ("s2".to_string(), 3)]
            .into_iter()
            .collect();
        assert_eq!(res, exp);
    }

    // #[test]
    pub fn testing2() {
        let t2g = "/home/michi/bus_testing/transcripts_to_genes.txt";
        // let b1 = BusFolder::new("/home/michi/bus_testing/bus_output/", t2g);
        let b1 = BusFolder::new("/home/michi/bus_testing/bus_output_short/");
        let b2 = BusFolder::new("/home/michi/bus_testing/bus_output_short/");

        // let hashmap = HashMap::from([
        //     ("full".to_string(), "/home/michi/bus_testing/bus_output/output.corrected.sort.bus".to_string()),
        //     ("short".to_string(), "/home/michi/bus_testing/bus_output_short/output.corrected.sort.bus".to_string())
        // ]);
        let hashmap = HashMap::from([("full".to_string(), b1), ("short".to_string(), b2)]);

        let s = make_fingerprint_histogram(hashmap, t2g);

        s.to_csv("/tmp/testing2.csv");
        println!("{:?}", s)
    }



    #[test]
    pub fn test_csv_read_write() {
        let order = vec!["s1", "s2"]
            .into_iter()
            .map(|x| x.to_string())
            .collect();
        let mut histogram: HashMap<Vec<u32>, usize> = HashMap::new();

        histogram.insert(vec![1, 0], 1);
        histogram.insert(vec![0, 1], 1);
        histogram.insert(vec![1, 1], 1);

        let fph = FingerprintHistogram {
            samplenames: order,
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
    pub fn test_estimate_sihr() {
        let order = vec!["s1", "s2"]
            .into_iter()
            .map(|x| x.to_string())
            .collect();
        let mut histogram: HashMap<Vec<u32>, usize> = HashMap::new();

        histogram.insert(vec![1, 0], 1);
        histogram.insert(vec![0, 1], 1);
        histogram.insert(vec![1, 1], 1);

        let fph = FingerprintHistogram {
            samplenames: order,
            histogram,
        };

        let p = fph.estimate_sihr();
        println!("{}", p);

        assert_almost_eq!(p, 0.5, 0.001);
    }

    #[test]
    fn test_sihr_real(){
        let fph = FingerprintHistogram::from_csv("/home/michi/Dropbox/rustphantompurger/IR56_57_phantom.csv");
        let sihr = fph.estimate_sihr();
        // check that the SIHR is still the same
        insta::assert_yaml_snapshot!(sihr, @r###"
        ---
        0.0016278754068794856
        "###);
    }

    #[test]
    fn test_zr_real(){
        let fph = FingerprintHistogram::from_csv("/home/michi/Dropbox/rustphantompurger/IR56_57_phantom.csv");
        // nubmer of non-chimeric reads
        let z_r = fph.get_z_r();
        
        // need to enforce hashmap sorting! otherwise the resulting Hashmap (Even if same)
        // will differ from the snapshot on disk
        let mut settings = insta::Settings::clone_current();
        settings.set_sort_maps(true);
        settings.bind(|| {
            // runs the assertion with the changed settings enabled
            insta::assert_yaml_snapshot!(z_r);
        });
    }

    #[test]
    fn test_mr_real(){
        let fph = FingerprintHistogram::from_csv("/home/michi/Dropbox/rustphantompurger/IR56_57_phantom.csv");
        // nubmer of total reads
        let m_r = fph.get_m_r();
        
        // need to enforce hashmap sorting! otherwise the resulting Hashmap (Even if same)
        // will differ from the snapshot on disk
        let mut settings = insta::Settings::clone_current();
        settings.set_sort_maps(true);
        settings.bind(|| {
            // runs the assertion with the changed settings enabled
            insta::assert_yaml_snapshot!(m_r);
        });
    }
}

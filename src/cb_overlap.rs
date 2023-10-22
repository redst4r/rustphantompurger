use std::{collections::HashMap, fs::File};
use flate2::{write::GzEncoder, Compression};
use itertools::Itertools;
use bustools::{io::BusReader, bus_multi::{CellIteratorMulti, CellUmiIteratorMulti}, utils::get_progressbar, iterators::{CbUmiGroupIterator, CellGroupIterator}};
use std::{
    hash::Hash,
    io::Write,
};
use crate::utils::valmap_ref;

#[derive(Eq, PartialEq, Hash, Ord, PartialOrd, Clone, Debug, Copy)]
pub struct CB(u64);

/// For each CB calcualte the number of UMIs per experiment/busfile.
/// If extensive swapping has occured, we'll see shared CBs across experiments
/// with correlated #umis
pub fn detect_cell_overlap(busfolders: &HashMap<String, String>, outfile: &str) {

    if !outfile.ends_with(".csv") && !outfile.ends_with(".csv.gz") {
        panic!("unknwon file extension. must be either .csv or .csv.gz")
    }

    // figure out size of iterators, just for progress bar!
    let cbs_per_file = valmap_ref(
        |busfile| {
            println!("determine size of iterator {busfile}");
            BusReader::new(busfile).groupby_cb().count()
        },
        busfolders,
    );

    println!("total records {:?}", cbs_per_file);
    let total: usize = cbs_per_file.values().sum();

    let samplenames: Vec<String> = busfolders.keys().cloned().collect();
    let multi_iter = CellIteratorMulti::new(busfolders);
    let mut result: HashMap<CB, Vec<usize>> = HashMap::new();

    let bar = get_progressbar(total as u64);

    for (i, (c, record_dict)) in multi_iter.enumerate() {
        let mut entry: Vec<usize> = Vec::new();
        for s in samplenames.iter() {
            let numi = match record_dict.get(s) {
                Some(records) => records.iter().map(|r| r.UMI).unique().count(),
                None => 0,
            };
            entry.push(numi)
        }
        result.insert(CB(c), entry);

        if i % 10_000 == 0 {
            // cells iterations are usually rather small, i.e. millions, update more reg
            bar.inc(10_000);
        }
    }

    if outfile.ends_with(".csv") {
        // write to file
        // TODO could be inlined into the above code to instantly write
        let mut fh = File::create(outfile).unwrap();
        let mut header = samplenames.join(",");
        header.push_str(",CB");
        writeln!(fh, "{}", header).unwrap();

        for (cid, numis) in result.iter() {
            // concat with commas
            let mut s = numis
                .iter()
                .map(|i| i.to_string())
                .collect::<Vec<String>>()
                .join(",");
            s.push_str(&format!(",{}", cid.0));
            writeln!(fh, "{}", s).unwrap();
        }
    } else if outfile.ends_with(".csv.gz") {
        // write compressed
        let fh = File::create(&format!("{outfile}.gz")).unwrap();
        let mut e = GzEncoder::new(fh, Compression::default());
        let mut header = samplenames.join(",");
        header.push_str(",CB\n");
        e.write_all(header.as_bytes()).unwrap();
        for (cid, numis) in result.iter() {
            // concat with commas
            let mut s = numis
                .iter()
                .map(|i| i.to_string())
                .collect::<Vec<String>>()
                .join(",");
            s.push_str(&format!(",{}\n", cid.0));
            e.write_all(s.as_bytes()).unwrap();
        }
    } else {
        panic!("unknwon file extension. must be either .csv or .csv.gz")
    }
}


pub fn detect_overlap(busfolders: HashMap<String, String>) -> HashMap<Vec<String>, usize> {
    // deprecated
    // mesures the number of ovelapping CUGs across the experiments

    let mut total = 0;
    // TODO: this doesnt check if the EC overlaps
    for v in busfolders.values() {
        println!("determine size of iterator");
        let total_records = BusReader::new(v).groupby_cbumi().count();
        if total < total_records {
            total = total_records
        }
    }
    println!("total records {}", total);

    let multi_iter = CellUmiIteratorMulti::new(&busfolders);
    let bar = get_progressbar(total as u64);
    let mut counter: HashMap<Vec<String>, usize> = HashMap::new();

    for (i, ((_cb, _umi), record_dict)) in multi_iter.enumerate() {
        let mut the_set: Vec<String> = record_dict.keys().cloned().collect();
        the_set.sort();
        let val = counter.entry(the_set).or_insert(0);
        *val += 1;

        if i % 1000000 == 0 {
            bar.inc(1000000);
        }
    }
    counter
}


#[cfg(test)]
mod test {
    use bustools::io::{setup_busfile, BusRecord};
    use super::*;
    
    #[test]
    fn test_detect_overlap() {
        let r1 =BusRecord{CB: 0, UMI: 1, EC: 0, COUNT: 2, FLAG: 0};
        let r2 = BusRecord{CB: 0, UMI: 1, EC: 0, COUNT: 3, FLAG: 0};
        let (fname1, _d1) = setup_busfile(&vec![r1]);
        let (fname2, _d2) = setup_busfile(&vec![r2]);
        let busfolders = HashMap::from([
            ("s1".to_string(), fname1.clone()),
            ("s2".to_string(), fname2.clone()),
        ]);
        println!("{}", fname1);
        detect_cell_overlap(&busfolders, "/tmp/overlap.csv");
        detect_cell_overlap(&busfolders, "/tmp/overlap.csv.gz");
    }

    // #[test]
    pub fn test_detect_cell_overlap() {
        // let t2g = "/home/michi/bus_testing/transcripts_to_genes.txt";
        // let b1 = BusFolder::new("/home/michi/bus_testing/bus_output/", t2g);
        // let b1 = BusFolder::new("/home/michi/bus_testing/bus_output_short/", t2g);
        // let b2 = BusFolder::new("/home/michi/bus_testing/bus_output_short/", t2g);
        // let busfolders = HashMap::from([
        //     ("full".to_string(), b1),
        //     ("short".to_string(), b2)
        // ]);

        let busfolders = HashMap::from([
            (
                "full".to_string(),
                "/home/michi/bus_testing/bus_output/output.corrected.sort.bus".to_string(),
            ),
            (
                "short".to_string(),
                "/home/michi/bus_testing/bus_output_short/output.corrected.sort.bus".to_string(),
            ),
        ]);

        detect_cell_overlap(&busfolders, "/tmp/test_detect_cell_overlap.csv")
    }
}
use std::collections::{HashMap, HashSet};
use rustbustools::{io::BusRecord, disjoint::DisjointSubsets, consistent_genes::Ec2GeneMapper};
use crate::phantompurger::{Fingerprint, groupby_gene_across_samples, make_fingerprint_simple};

/*
some old code to assemble fingerprints,
deprecated
 */
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


// #[inline(never)]
fn groupby_gene_even_simpler(record_dict: HashMap<String,BusRecord>, ecmapper_dict: &HashMap<String, &Ec2GeneMapper>) -> Vec<HashMap<String, BusRecord>> {
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
    use std::collections::HashMap;

    use rustbustools::{io::BusRecord, consistent_genes::Ec2GeneMapper};

    use crate::phantompurger::create_dummy_ec;

    use super::{groupby_gene_even_simpler};

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

}
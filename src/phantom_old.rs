fn groupby_gene(record_dict: HashMap<String,Vec<BusRecord>>, ecmapper_dict: &HashMap<String, Ec2GeneMapper>) 
    -> Vec<HashMap<String, BusRecord>>
    {
    // multiple experiments yielded the same CB/UMI
    // turn this into a list of record_dicts, where the gene is consistent

    let mut records:Vec<BusTmp> = Vec::new();
    for (sname, recordlist) in record_dict{
        for r in recordlist{
            records.push(
                BusTmp { cb: r.CB, umi: r.UMI, ec: r.EC, count: r.COUNT, flag: r.FLAG, samplename: sname.clone() }
            )
        }
    }
    // if records.len() == 1{
    //     return record_dict
    //     // return
    // };

    // get genes for each record.
    let mut genes: HashMap<BusTmp, &HashSet<u32>> = HashMap::new();
    for b in records{
        let ecmapper = ecmapper_dict.get(&b.samplename).unwrap();
        let g = ecmapper.get_genes(b.ec);
        genes.insert(b, g);
    }

    // also build the disjoint set based on the genes
    let mut disjoint_set = DisjointSubsets::new();
    for (b, gset) in genes.drain(){
        disjoint_set.add(b.to_string(), gset.clone());
    }

    let mut emit_vector: Vec<HashMap<String, BusRecord>> = Vec::new();
    for bus_string_concat in disjoint_set.disjoint_sets.keys(){
        // these are now all records across the samples that map to the same gene
        let bus_tuple: Vec<BusTmp> = bus_string_concat.split('_').map(|bstring|BusTmp::parse(bstring).unwrap()).collect();

        let mut emited_dict: HashMap<String, BusRecord> = HashMap::new();
        for b in bus_tuple{

            // if we already have a record in that sample, just update the count
            if emited_dict.contains_key(&b.samplename){
                let brecord = emited_dict.get_mut(&b.samplename).unwrap();
                brecord.COUNT += b.count;
            } 
            else{
                let brecord = BusRecord{ CB: b.cb, UMI: b.umi, EC:0, COUNT: b.count, FLAG: b.flag};
                emited_dict.insert(b.samplename, brecord);
            }
        }
        emit_vector.push(emited_dict);
    }
    emit_vector
}

// #[inline(never)]
pub fn groupby_gene_simple(record_dict: HashMap<String,BusRecord>, ecmapper_dict: &HashMap<String, &Ec2GeneMapper>) -> Vec<HashMap<String, BusRecord>>
{
    // assuming no multimapped reads, hence record dict contains only single entries (not a list)
    // this uses the awkward conversion to BusTmp


    // check if the genes are consistent across samples
    // if so yield the record as is
    // otherwise split into multiple records
    if record_dict.len() == 1{
        return vec![record_dict];
    }
    
    let mut records:Vec<BusTmp> = Vec::new();
    for (sname, r) in record_dict{
        records.push(
                BusTmp { cb: r.CB, umi: r.UMI, ec: r.EC, count: r.COUNT, flag: r.FLAG, samplename: sname.clone() }
            )
    }

    let mut genes: HashMap<BusTmp, &HashSet<u32>> = HashMap::new();
    for b in records{
        let ecmapper = ecmapper_dict.get(&b.samplename).unwrap();
        let g = ecmapper.get_genes(b.ec);
        genes.insert(b, g);
    }

    let mut disjoint_set = DisjointSubsets::new();
    for (b, gset) in genes.drain(){
        disjoint_set.add(b.to_string(), gset.clone());
    }

    let mut emit_vector: Vec<HashMap<String, BusRecord>> = Vec::new();
    for bus_string_concat in disjoint_set.disjoint_sets.keys(){
        // these are now all records across the samples that map to the same gene
        let bus_tuple: Vec<BusTmp> = bus_string_concat.split(SEPARATOR).map(|bstring|BusTmp::parse(bstring).unwrap()).collect();

        let mut emited_dict: HashMap<String, BusRecord> = HashMap::new();
        for b in bus_tuple{

            if emited_dict.contains_key(&b.samplename){
                panic!("cant happen, each sample only has one record")
            } 
            else{
                let brecord = BusRecord{ CB: b.cb, UMI: b.umi, EC:0, COUNT: b.count, FLAG: b.flag};
                emited_dict.insert(b.samplename, brecord);
            }
        }
        emit_vector.push(emited_dict);
    }
    emit_vector
}

#[derive(Eq, PartialEq, Debug, Hash)]
struct BusTmp {cb: u64, umi: u64, ec: u32, count:u32, flag:u32, samplename:String }
impl BusTmp {
    fn parse(item: &str) -> Option<Self> {

        let s: Vec<_> = item.split("@@").collect();
        if s.len() == 6{
            let cb:u64 = s[0].parse().unwrap();
            let umi:u64 = s[1].parse().unwrap();
            let ec: u32 = s[2].parse().unwrap();
            let count: u32 = s[3].parse().unwrap();
            let flag:u32 = s[4].parse().unwrap();
            let samplename = s[5].to_string();
            Some(BusTmp {cb , umi , ec, count, flag, samplename})
        }
        else{
            None
        }
    }

    fn to_string(&self) -> String{
        format!("{}@@{}@@{}@@{}@@{}@@{}", self.cb, self.umi, self.ec, self.count, self.flag, self.samplename)    }
}

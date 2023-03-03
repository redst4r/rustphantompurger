use std::{collections::{HashMap, HashSet}, hash::Hash};

pub struct DisjointSubsets<T>{
    pub disjoint_sets: HashMap<String, HashSet<T>>
}

pub const SEPARATOR: &str ="@SEP@";

impl <T:Hash+Eq> DisjointSubsets<T>{
    pub fn new() -> Self{
        DisjointSubsets{
            disjoint_sets: HashMap::new()
        }
    }

    pub fn add(&mut self, setname:String, aset:HashSet<T>){

        if self.disjoint_sets.contains_key(&setname){
            panic!("inserting into existing setname")
        }
        // check for any existing set that shares a member with aset
        let mut candidates: Vec<String> = Vec::new();
        for (name, the_set) in &self.disjoint_sets{
            if set_overlap(&aset, the_set){
                candidates.push(name.clone());
            }
        }

        // aset links everything in the candidates into a single cell
        // create that new set!
        let mut newset:HashSet<T> = HashSet::new();
        let mut newname = String::new();
        for c in candidates{
            let mut s = self.disjoint_sets.remove(&c).unwrap(); //poing the old set out
            for el in s.drain(){
                newset.insert(el);
            }
            newname.push_str(&c);
            // newname.push('_');
            newname.push_str(SEPARATOR);
        }

        // add the set itself
        for el in aset{
            newset.insert(el);
        }
        newname.push_str(&setname);

        // insert the newly made set
        self.disjoint_sets.insert(newname, newset);
    }

    pub fn get_disjoint_set_ids(&self) -> Vec<Vec<String>>{
        // return the list of disjoint sets, for each set, report it by its element's IDs
        let mut setlist: Vec<Vec<String>> = Vec::with_capacity(self.disjoint_sets.len());
        for id_str in self.disjoint_sets.keys(){
            let s: Vec<String> = id_str.split(SEPARATOR).map(|x|x.to_string()).collect();
            setlist.push(s);
        }
        setlist
    }

    pub fn get_disjoint_set_ids_and_set(&self) -> Vec<(Vec<String>, &HashSet<T>)>{
        // return the list of disjoint sets, for each set, report it by its element's IDs and the set
        let mut setlist: Vec<(Vec<String>, &HashSet<T>)> = Vec::with_capacity(self.disjoint_sets.len());
        for id_str in self.disjoint_sets.keys(){
            let s: Vec<String> = id_str.split(SEPARATOR).map(|x|x.to_string()).collect();
            let theset = self.disjoint_sets.get(id_str).unwrap();
            setlist.push((s, theset));
        }
        setlist
    }

}

fn set_overlap<T: Hash+Eq>(aset: &HashSet<T>, bset: &HashSet<T>) -> bool{
    for el in aset{
        if bset.contains(el){
            return true
        }
    }
    false
}

#[cfg(test)]
mod testing{
    use std::collections::{HashSet, HashMap};

    use crate::disjoint::DisjointSubsets;

    #[test]
    fn testing_disjoint(){

        let a: HashSet<_> = vec!["A", "B","C"].into_iter().collect();
        let b: HashSet<_> = vec!["D", "E"].into_iter().collect();
        let c: HashSet<_> = vec!["B", "D"].into_iter().collect();
        let d: HashSet<_> = vec!["Z"].into_iter().collect();

        let mut ds = DisjointSubsets::new();
        ds.add("set1".to_string(), a);
        ds.add("set2".to_string(), b);


        let mut expected: HashMap<String, HashSet<&str>> = HashMap::new();
        let set1 = vec!["A", "B", "C"].into_iter().collect::<HashSet<&str>>();
        expected.insert("set1".to_string(), set1);
        let set2 = vec!["D", "E"].into_iter().collect::<HashSet<&str>>();
        expected.insert("set2".to_string(), set2);

        assert_eq!(
            ds.disjoint_sets,
            expected
        );

        // println!("{:?}", ds.disjoint_sets);

        ds.add("set3".to_string(), c);

        let mut expected: HashMap<String, HashSet<&str>> = HashMap::new();
        let set1 = vec!["A", "B", "C", "D","E"].into_iter().collect::<HashSet<&str>>();
        expected.insert("set1@SEP@set2@SEP@set3".to_string(), set1);
        assert_eq!(
            ds.disjoint_sets,
            expected
        );
        // println!("{:?}", ds.disjoint_sets);


        ds.add("set4".to_string(), d);

        let mut expected: HashMap<String, HashSet<&str>> = HashMap::new();
        let set1 = vec!["A", "B", "C", "D","E"].into_iter().collect::<HashSet<&str>>();
        expected.insert("set1@SEP@set2@SEP@set3".to_string(), set1);
        let set2 = vec!["Z"].into_iter().collect::<HashSet<&str>>();
        expected.insert("set4".to_string(), set2);    
        // println!("{:?}", ds.disjoint_sets);
        assert_eq!(
            ds.disjoint_sets,
            expected
        );
    }
}
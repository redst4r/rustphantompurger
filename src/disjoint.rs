use std::{
    collections::HashSet,
    hash::Hash,
};
use std::collections::HashMap;

fn set_overlap<T: Hash + Eq>(aset: &HashSet<T>, bset: &HashSet<T>) -> bool {
    for el in aset {
        if bset.contains(el) {
            return true;
        }
    }
    false
}

#[derive(Debug)]
pub struct DisjointSubsets<T> {
    pub disjoint_sets: HashMap<String, HashSet<T>>,
}

pub const SEPARATOR: &str = "@SEP@";

impl<T: Hash + Eq> DisjointSubsets<T> {
    pub fn new() -> Self {
        DisjointSubsets {
            disjoint_sets: HashMap::new(),
        }
    }

    pub fn add(&mut self, setname: String, aset: HashSet<T>) {
        if self.disjoint_sets.contains_key(&setname) {
            panic!("inserting into existing setname")
        }
        // check for any existing set that shares a member with aset
        let mut candidates: Vec<String> = Vec::new();
        for (name, the_set) in &self.disjoint_sets {
            if set_overlap(&aset, the_set) {
                candidates.push(name.clone());
            }
        }

        // aset links everything in the candidates into a single cell
        // create that new set!
        let mut newset: HashSet<T> = HashSet::new();
        let mut newname = String::new();
        for c in candidates {
            let mut s = self.disjoint_sets.remove(&c).unwrap(); //poing the old set out
            for el in s.drain() {
                newset.insert(el);
            }
            newname.push_str(&c);
            // newname.push('_');
            newname.push_str(SEPARATOR);
        }

        // add the set itself
        for el in aset {
            newset.insert(el);
        }
        newname.push_str(&setname);

        // insert the newly made set
        self.disjoint_sets.insert(newname, newset);
    }

    pub fn get_disjoint_set_ids(&self) -> Vec<Vec<String>> {
        // return the list of disjoint sets, for each set, report it by its element's IDs
        let mut setlist: Vec<Vec<String>> = Vec::with_capacity(self.disjoint_sets.len());
        for id_str in self.disjoint_sets.keys() {
            let s: Vec<String> = id_str.split(SEPARATOR).map(|x| x.to_string()).collect();
            setlist.push(s);
        }
        setlist
    }

    pub fn get_disjoint_set_ids_and_set(&self) -> Vec<(Vec<String>, &HashSet<T>)> {
        // return the list of disjoint sets, for each set, report it by its element's IDs and the set
        let mut setlist: Vec<(Vec<String>, &HashSet<T>)> =
            Vec::with_capacity(self.disjoint_sets.len());
            
        for id_str in self.disjoint_sets.keys() {
            let s: Vec<String> = id_str.split(SEPARATOR).map(|x| x.to_string()).collect();
            let theset = self.disjoint_sets.get(id_str).unwrap();
            setlist.push((s, theset));
        }
        setlist
    }
}



#[cfg(test)]
mod test {
    use super::*;
    fn vec2set<T: Eq + Hash>(x: Vec<T>) -> HashSet<T> {
        x.into_iter().collect::<HashSet<T>>()
    }

    // use super::vec2set;

    use crate::disjoint::DisjointSubsets;
    use std::collections::{HashMap, HashSet};

    #[test]
    fn test_get_disjoint_set_ids() {
        let a = vec2set(vec!["A", "B", "C"]);
        let b = vec2set(vec!["D", "E"]);
        let c = vec2set(vec!["B", "D"]);

        let mut ds = DisjointSubsets::new();
        ds.add("set1".to_string(), a);
        ds.add("set2".to_string(), b);

        let disjoint_setids = ds.get_disjoint_set_ids();
        /*
        at this point it looks liek
        setnames:  {{"set1"}, {"set2"}}
        setvalues: {{A, B, C}, {D, E }}
         */

        // a bit tricky with the nested sets!!
        // actually rust doesnt allow this, HashSet doesnt implement Hash
        /*let exp = 
        vec2set(vec![
            vec2set(vec!["set1".to_string()]), 
            vec2set(vec!["set2".to_string()])
        ]);
        let obs = 
        vec2set(vec![
            vec2set(disjoint_setids[0].clone()), 
            vec2set(disjoint_setids[1].clone())
        ]);
        assert_eq!(exp, obs);*/
        assert_eq!(disjoint_setids.len(), 2);       


        // add another set, linking the two
        /*
        at this point it looks liek
        setnames:  {{"set1", "set2", "set3"}}
        setvalues: {{A, B, C D, E }}
         */
        ds.add("set3".to_string(), c);
        println!("{:?}", ds);
        assert_eq!(disjoint_setids.len(), 3);       
        assert_eq!(
            vec2set(disjoint_setids[0].clone()),
            vec2set(vec!["set1".to_string(), "set2".to_string(), "set3".to_string()])
        );
    }

    #[test]
    fn testing_disjoint() {
        let a = vec2set(vec!["A", "B", "C"]);
        let b = vec2set(vec!["D", "E"]);
        let c = vec2set(vec!["B", "D"]);
        let d = vec2set(vec!["Z"]);

        let mut ds = DisjointSubsets::new();
        ds.add("set1".to_string(), a);
        ds.add("set2".to_string(), b);

        let mut expected: HashMap<String, HashSet<&str>> = HashMap::new();
        // let set1 = vec2set(vec!["vec2set(A", "B", "C"]);
        let set1 = vec2set(vec!["A", "B", "C"]);
        expected.insert("set1".to_string(), set1);
        let set2 = vec2set(vec!["D", "E"]);
        expected.insert("set2".to_string(), set2);

        assert_eq!(ds.disjoint_sets, expected);

        // println!("{:?}", ds.disjoint_sets);

        // adding set3 will connect everything!
        ds.add("set3".to_string(), c);
        assert_eq!(ds.disjoint_sets.len(), 1);

        // this test doesnt work as the ID "set1@SEP@set2@SEP@set3" is randomly order
        // let mut expected: HashMap<String, HashSet<&str>> = HashMap::new();
        // let set1 = vec2set(vec!["A", "B", "C", "D", "E"]);
        // expected.insert("set1@SEP@set2@SEP@set3".to_string(), set1);
        // assert_eq!(ds.disjoint_sets, expected);

        // added set 4: disjoint from everything
        ds.add("set4".to_string(), d);

        let mut expected: HashMap<String, HashSet<&str>> = HashMap::new();
        let set1 = vec2set(vec!["A", "B", "C", "D", "E"]);
        expected.insert("set1@SEP@set2@SEP@set3".to_string(), set1);
        let set2 = vec2set(vec!["Z"]);
        expected.insert("set4".to_string(), set2);
        // println!("{:?}", ds.disjoint_sets);
        assert_eq!(ds.disjoint_sets, expected);
    }
}
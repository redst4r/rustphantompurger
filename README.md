# Phantompurger in Rust

A Rust version of the [PhantompPurger algorithm](https://www.nature.com/articles/s41467-020-16522-z) to remove chimeric molecules from scRNAseq data. 

**Far from complete**, in particular we're not implementing any of the statistical approaches to decide which molecules are chimeras.
Rather, this crate (quickly!) creates the input tables to feed into the PhantomPurger statistical framework from `bus`-files (i.e. kallisto output).


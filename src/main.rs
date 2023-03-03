use std::collections::HashMap;
use rustbustools::io::BusFolder;
use rustphantompurger::phantompurger;
use clap::{self, Parser, Subcommand, Args};

#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
struct Cli {  
    /// Path to output file
    #[clap(short ='o', long = "output")] 
    output: String,    

    #[clap(subcommand)]
    command: MyCommand
}

#[allow(non_camel_case_types)]
#[derive(Subcommand)]
enum MyCommand {
    phantomFingerprint(PhantomArgs),
    phantomCB(PhantomCBArgs),
    phantomEstimate(PhantomEstArgs),
    phantomFilter(PhantomFilterArgs),

}
#[derive(Args)]
struct PhantomArgs{
    #[clap(long= "infolders")]
    busfolders: String,
    #[clap(long= "t2g")] 
    t2g: String, 
}

#[derive(Args)]
struct PhantomCBArgs{
    #[clap(long= "infolders")]
    busfolders: String,
}

#[derive(Args)]
struct PhantomEstArgs{
    #[clap(long= "csv")]
    phantomcsv: String,
}

#[derive(Args)]
struct PhantomFilterArgs{
    #[clap(long= "csv")]
    phantomcsv: String,

    #[clap(long= "infolders")]
    infolders: String,

    #[clap(long= "outfiles")]
    outfiles: String,

    #[clap(long= "removed")]
    removed: String,    
    
    #[clap(long= "threshold")]
    threshold: f64,

    #[clap(long= "t2g")] 
    t2g: String, 
}

fn main() {
    let cli = Cli::parse();
    match cli.command{
        MyCommand::phantomFingerprint(args) => {
            println!("Doing phantom");
            let named_infolders= parse_key_value_args(&args.busfolders);
            println!("Infiles: {:?}",named_infolders);

            let busfolder_dict = named_infolders.into_iter()
                .map(|(sname, folder)|
                    (sname, BusFolder::new(&folder, &args.t2g)))
                .collect();

            let histo = phantompurger::make_fingerprint_histogram(busfolder_dict);
            histo.to_csv(&cli.output);
        }
        MyCommand::phantomCB(args) => {
            println!("Doing phnatom CB overlap");
            let named_infolders= parse_key_value_args(&args.busfolders);
            phantompurger::detect_cell_overlap(named_infolders, &cli.output);
        }
        MyCommand::phantomEstimate(args) => {
            println!("Doing phantom SIHR estimation");

            let fph = phantompurger::FingerprintHistogram::from_csv(&args.phantomcsv);

            let p = fph.estimate_sihr();
            println!("Estimated SIHR as {}", p);
        }
        MyCommand::phantomFilter(args) => {

            let named_infolders= parse_key_value_args(&args.infolders);
            let named_outfiles= parse_key_value_args(&args.outfiles);
            let named_removed= parse_key_value_args(&args.removed);
            println!("Infiles: {:?}",named_infolders);
            println!("outfiles: {:?}",named_outfiles);
            println!("named_removed: {:?}",named_removed);

            let fph = phantompurger::FingerprintHistogram::from_csv(&args.phantomcsv);

            println!("Building posterior");
            let posterior = phantompurger::PhantomPosterior::new(&fph);

            println!("Building busfolder dicts");
            let inputfolder_dict = named_infolders.iter()
                .map(|(name, folder)|
                    (name.clone(),  BusFolder::new(folder, &args.t2g)))
                .collect();

            println!("Filtering");

            posterior.filter_busfiles(
                inputfolder_dict, 
                named_outfiles, 
                named_removed,
                args.threshold
            )
        }        
    }
}

fn parse_key_value_args(s: &str) -> HashMap<String,String>{

    // we supply key:value args on the command line via this syntax:
    // k1:v1 k2:v2 ....
    // parse into dict
    let mut the_map = HashMap::new();

    for pair in  s.split(' '){
        let kv:Vec<&str> = pair.split(':').collect();
        assert!(kv.len()==2, "{pair}");
        let k = kv[0].to_string();
        let v = kv[1].to_string();

        the_map.insert(k, v);
    }
    the_map

}

/*
cargo run --release -- --output /dev/null phantom-estimate --csv /home/michi/Dropbox/rustphantompurger/IR56_57_phantom.csv

 */
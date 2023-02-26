use itertools::izip;
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
    #[clap()]
    busfolders: Vec<String>,
    #[clap(long= "t2g")] 
    t2g: String, 
}

#[derive(Args)]
struct PhantomCBArgs{
    #[clap()]
    busfolders: Vec<String>,
}

#[derive(Args)]
struct PhantomEstArgs{
    #[clap()]
    phantomcsv: String,
}

#[derive(Args)]
struct PhantomFilterArgs{
    #[clap()]
    phantomcsv: String,
    #[clap()]
    infolders: Vec<String>,
    #[clap()]
    outfiles: Vec<String>,
    #[clap(long= "t2g")] 
    t2g: String, 

}

fn main() {
    let cli = Cli::parse();
    match cli.command{
        MyCommand::phantomFingerprint(args) => {
            println!("Doing phantom");

            let busfolder_dict = args.busfolders.into_iter()
                .map(|b|(b.clone(), BusFolder::new(&b, &args.t2g)))
                .collect();

            let histo = phantompurger::make_fingerprint_histogram(busfolder_dict);
            histo.to_csv(&cli.output);
        }
        MyCommand::phantomCB(args) => {
            println!("Doing phnatom CB overlap");

            let busfolder_dict = args.busfolders.into_iter()
                .map(|b|(b.clone(), b))
                .collect();
            phantompurger::detect_cell_overlap(busfolder_dict, &cli.output);
        }
        MyCommand::phantomEstimate(args) => {
            println!("Doing phantom SIHR estimation");

            let fph = phantompurger::FingerprintHistogram::from_csv(&args.phantomcsv);

            let p = fph.estimate_sihr();
            println!("Estimated SIHR as {}", p);
        }
        MyCommand::phantomFilter(args) => {

            let fph = phantompurger::FingerprintHistogram::from_csv(&args.phantomcsv);
            let posterior = phantompurger::PhantomPosterior::new(&fph);
            let inputfolder_dict = args.infolders.iter()
                .map(|b|(b.clone(),  BusFolder::new(&b, &args.t2g)))
                .collect();

            let output_busfolders = izip!(
                args.infolders.iter(), 
                args.outfiles.iter()
            )
            .map(|(name, outfile)|(name.clone(), outfile.clone()))
            .collect();

            posterior.filter_busfiles(inputfolder_dict, output_busfolders)
        }        
    }
}
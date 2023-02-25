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
    phantom(PhantomArgs),
    phantomCB(PhantomCBArgs),
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

fn main() {
    let cli = Cli::parse();
    match cli.command{
        MyCommand::phantom(args) => {
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
    }
}
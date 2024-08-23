use clap::Parser;
use rust_htslib::bam::{Read, Reader, Record, Header};
use std::fs::File;
use std::io::{self, Write};
use bam_aligner::bam_subset::BamSubset;
use num_cpus;

/// Simple program to assemble a contig sequence from a BAM file
#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Cli {
    /// BAM file path
    #[arg(short, long)]
    bam: String,

    /// Output FASTA file name
    #[arg(short, long)]
    outfile: String,

    /// number of threads to use
    #[arg(short, long)]
    num_threads: Option<usize>,
}



fn main() -> io::Result<()> {
    // Parse command line arguments
    let args = Cli::parse();

    let num_threads = match args.num_threads{
        Some(n) => {
            if n < 1{
                1
            } 
            else {
                n
            }
        },
        None => num_cpus::get(),
    };

    // Read BAM file and collect reads
    let records = read_bam_file(&args.bam, num_threads);

    println!("I collected info on {} cells", records.len() );

    // Assemble contig from reads
    let contig_sequences = records.assemble_contigs();

    // Write the contig sequence to a FASTA file
    write_fasta(&args.outfile, &contig_sequences)?;

    println!("Contig assembled and written to {}", args.outfile);

    Ok(())
}


fn read_bam_file(file_path: &str, num_threads:usize ) -> BamSubset {
    // 1. Open the BAM file using rust-htslib's Reader
    let mut bam = Reader::from_path(file_path).expect("Error opening BAM file");

    // 2. Initialize an empty vector to store BAM records
    let mut records = BamSubset::new(num_threads);

    // 3. Iterate over each read in the BAM file
    for result in bam.records() {

        // 4. Handle potential errors when reading a record
        let record = result.expect("Failed to read BAM record");
        let tid = record.tid();
        let this_start = record.pos() as u32;
        let this_end = record.cigar().end_pos() as u32;

        records.push( record.clone() );
        
    }

    // 7. Return the BamSubset for further processing
    records
}



// Function to write the assembled contig to a FASTA file
fn write_fasta(file_path: &str, contig_sequence: &Vec< String>) -> io::Result<()> {
    let mut file = File::create(file_path)?;
    for fasta in contig_sequence{
        writeln!(file, "{fasta}", )?;
    }
    Ok(())
}

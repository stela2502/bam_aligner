use clap::Parser;
use rust_htslib::bam::{Read, Reader, Record, Header};
use std::fs::File;
use std::io::{self, Write};


/// Simple program to assemble a contig sequence from a BAM file
#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Cli {
    /// BAM file path
    #[arg(short, long)]
    bam: String,

    /// Output FASTA file name
    #[arg(short, long)]
    contig: String,

    /// The fasta accession id
    #[arg(short, long)]
    accession: String,
}



fn main() -> io::Result<()> {
    // Parse command line arguments
    let args = Cli::parse();

    // Read BAM file and collect reads
    let records = read_bam_file(&args.bam);

    // Assemble contig from reads
    let contig_sequence = assemble_contig(records);

    // Write the contig sequence to a FASTA file
    write_fasta(&args.contig, &args.accession, &contig_sequence)?;

    println!("Contig assembled and written to {}", args.contig);

    Ok(())
}


fn read_bam_file(file_path: &str) -> Vec<Record> {
    // 1. Open the BAM file using rust-htslib's Reader
    let mut bam = Reader::from_path(file_path).expect("Error opening BAM file");

    // 2. Initialize an empty vector to store BAM records
    let mut records = Vec::new();

    // 3. Iterate over each read in the BAM file
    for result in bam.records() {
        // 4. Handle potential errors when reading a record
        let record = result.expect("Error reading BAM record");

        // 5. Apply any necessary filters (e.g., minimum quality, specific regions)
        // Example: filter out unmapped reads
        if record.is_unmapped() {
            continue;
        }

        // 6. Push the valid record to the records vector
        records.push(record);
    }

    // 7. Return the vector of records for further processing
    records
}


fn assemble_contig(mut records: Vec<Record>) -> String {
    // 1. Sort the records by their start positions on the reference genome
    records.sort_by_key(|r| r.pos());

    // 2. Initialize an empty string to hold the consensus sequence
    let mut consensus_sequence = String::new();

    // 3. Initialize variables to track the current contig being assembled
    let mut current_contig_start: i64;
    let mut current_contig_end = 0;
    let mut contig_reads = Vec::new();

    // 4. Iterate through sorted records
    for record in records {
        let start_pos = record.pos();
        let end_pos = record.cigar().end_pos();

        // 5. Check if the read overlaps with the current contig
        if start_pos <= current_contig_end {
            // 6. If it overlaps, add it to the current contig's reads
            contig_reads.push(record);
            // Update the contig end position if this read extends beyond the current end
            if end_pos > current_contig_end {
                current_contig_end = end_pos;
            }
        } else {
            // 7. If it doesn't overlap, process the previous contig
            if !contig_reads.is_empty() {
                let consensus = build_consensus(&contig_reads);
                consensus_sequence.push_str(&consensus);
            }

            // 8. Start a new contig
            contig_reads.clear();
            contig_reads.push(record);
            current_contig_start = start_pos;
            current_contig_end = end_pos;
        }
    }

    // 9. Handle the last set of contig reads
    if !contig_reads.is_empty() {
        let consensus = build_consensus(&contig_reads);
        consensus_sequence.push_str(&consensus);
    }

    // 10. Return the assembled consensus sequence
    consensus_sequence
}

// Helper function to build consensus from a set of overlapping reads
fn build_consensus(contig_reads: &[Record]) -> String {
    // 1. Initialize a vector to hold consensus bases for each position
    let mut consensus_bases: Vec<char> = Vec::new();

    // 2. Determine the range of positions covered by the reads
    let start_pos = contig_reads[0].pos();
    let end_pos = contig_reads.iter().map(|r| r.cigar().end_pos()).max().unwrap();

    // 3. Iterate through each position within this range
    for pos in start_pos..end_pos {
        let mut base_counts = std::collections::HashMap::new();

        // 4. For each read that covers this position, count the occurrences of each base
        for read in contig_reads {

            if pos >= read.pos() && pos < read.cigar().end_pos() {
                // Convert position to index within the read's sequence
                let idx = (pos - read.pos()) as usize;
                // Access the base at this position
                let array = read.seq().as_bytes();
                if array.len() > idx{
                    let base = read.seq().as_bytes()[idx] as char;
                    *base_counts.entry(base).or_insert(0) += 1;
                }
                
            }
        }

        // 5. Determine the most common base at this position (consensus)
        if let Some((&consensus_base, _)) = base_counts.iter().max_by_key(|&(_, count)| count) {
            consensus_bases.push(consensus_base);
        }
    }

    // 6. Convert the vector of consensus bases to a string and return it
    consensus_bases.iter().collect()
}


// Function to write the assembled contig to a FASTA file
fn write_fasta(file_path: &str, acc:&str, contig_sequence: &str) -> io::Result<()> {
    let mut file = File::create(file_path)?;
    writeln!(file, ">{}", acc)?;
    writeln!(file, "{}", contig_sequence)?;

    Ok(())
}

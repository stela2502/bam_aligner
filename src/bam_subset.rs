use rust_htslib::bam::{Record};
use std::collections::HashSet;
use std::collections::HashMap;
use rust_htslib::bam::record::Aux;


pub struct BamSubset {
    records: HashMap<String, Vec::<Record>>,
    cells: HashSet<String>,
}

impl BamSubset{
    
    pub fn new() -> Self {
        Self{
            records: HashMap::new(),
            cells: HashSet::<String>::new(),
        }
    }

    pub fn push(&mut self, record:Record ) {

        if let Some(tag_value) = record.aux(b"CB").ok() {
            if let Aux::String(barcode) = tag_value {
                // Convert the barcode to a String and do something with it
                let record_barcode = barcode.to_string();
                // Now you can use record_barcode
                self.records.entry(record_barcode.clone())
                    .or_insert_with(Vec::new)
                    .push(record);
                self.cells.insert(record_barcode);
            }
        }else {
            eprintln!("I could not detect the barcode in this BAM record: {record:?}");
        }
    }

    pub fn len( &self) -> usize{
        self.cells.len()
    }

    pub fn assemble_contigs( &self ) -> Vec<(String, String)> {
        // 1. Sort the records by their start positions on the reference genome

        // 2. Initialize an empty fasta vector to hold the consensus sequences (Acc, Seq)
        let mut consensus_sequence = Vec::<(String, String)>::new();

        // 3. Initialize variables to track the current contig being assembled
        let mut current_contig_start: i64;
        let mut current_contig_end = 0;
        let mut contig_reads = Vec::new();
        let mut id: usize;
        // 4. Iterate through sorted records
        for (cell, records) in &self.records {
            id = 1;
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
                        let consensus = self.build_consensus(&contig_reads);
                        consensus_sequence.push( (format!("{cell}|transcript_{id}"), consensus) );
                        id +=1;
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
                let consensus = self.build_consensus(&contig_reads);
                consensus_sequence.push( (format!("{cell}|transcript_{id}"), consensus) );
            }
        }

        // 10. Return the assembled consensus sequence
        consensus_sequence
    }

    // Helper function to build consensus from a set of overlapping reads
    pub fn build_consensus( &self, contig_reads:&Vec<&Record> ) -> String {
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
}


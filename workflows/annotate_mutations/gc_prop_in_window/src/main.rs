use bio::io::fasta;
use bio::seq_analysis::gc::gc_content;
use clap::{Arg, Command};
use rust_htslib::bcf::{Read, Reader};
use std::cmp;
use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::str;

fn main() -> Result<(), Box<dyn Error>> {
    let matches = Command::new("VCF GC Calculator")
        .arg(
            Arg::new("vcf_path")
                .long("vcf")
                .value_name("PATH")
                .required(true)
                .help("path to a VCF file"),
        )
        .arg(
            Arg::new("fasta_path")
                .long("fasta")
                .value_name("PATH")
                .required(true)
                .help("path to the corresponding reference FASTA"),
        )
        .arg(
            Arg::new("window_size")
                .long("window")
                .value_name("SIZE")
                .required(true)
                .help("number of base pairs to construct a centered window around the position"),
        )
        .arg(
            Arg::new("output_path")
                .long("output")
                .value_name("PATH")
                .required(true)
                .help("path to an output TSV file to be created"),
        )
        .get_matches();

    let vcf_path = matches.get_one::<String>("vcf_path").unwrap();
    let fasta_path = matches.get_one::<String>("fasta_path").unwrap();
    let window_size: u64 = matches.get_one::<String>("window_size").unwrap().parse()?;
    let output_path = matches.get_one::<String>("output_path").unwrap();

    process_vcf(vcf_path, fasta_path, window_size, output_path)?;

    println!("GC proportions calculated and written to {}", output_path);
    Ok(())
}

/// Calculate the GC proportion in a window around all positions in a VCF file and output a TSV file
/// with the results.
///
/// # Arguments
///
/// * `vcf_path`: path to a VCF file
/// * `fasta_path`: path to the corresponding reference FASTA
/// * `window_size`: number of base pairs to construct a centered window around the position
/// * `output_path`: path to an output TSV file to be created
fn process_vcf(
    vcf_path: &str,
    fasta_path: &str,
    window_size: u64,
    output_path: &str,
) -> Result<(), Box<dyn Error>> {
    // open the VCF file for reading
    let mut vcf = Reader::from_path(vcf_path)?;

    // extract the header for contig name lookup
    let header = vcf.header().to_owned();

    // open the FASTA file for reading
    let fasta_path_buf = PathBuf::from(fasta_path);
    let mut fasta_handle = fasta::IndexedReader::from_file(&fasta_path_buf)?;

    // write output to a TSV file
    let output_file = File::create(output_path)?;
    let mut writer = BufWriter::new(output_file);
    writeln!(writer, "CHROM\tPOS\tREF\tALT\tGC_PROP")?;

    // iterate over each record in the VCF file
    for record_result in vcf.records() {
        let record = record_result?;

        // extract the chromosome, position, and reference sequence from the VCF record
        let chrom = str::from_utf8(header.rid2name(record.rid().unwrap()).unwrap()).unwrap();
        let pos = (1 + record.pos()) as u64; // 1-indexing the positions for consistency
        let ref_seq = str::from_utf8(record.alleles()[0])?;
        let alt_seq = str::from_utf8(record.alleles()[1])?;

        // determine the variant class (to potentially adjust the window center position later)
        let variant_class = if record.alleles()[1].len() < ref_seq.len() {
            "deletion"
        } else if record.alleles()[1].len() > ref_seq.len() {
            "insertion"
        } else {
            "substitution"
        };

        // calculate GC proportion
        let gc_proportion = calc_gc_proportion(
            chrom,
            pos,
            ref_seq,
            variant_class,
            window_size,
            &mut fasta_handle,
        )?;

        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{:.4}",
            chrom, pos, ref_seq, alt_seq, gc_proportion
        )?;
    }

    Ok(())
}

/// Calculate the GC proportion in a window centered around a genomic position.
///
/// # Arguments
///
/// * `chrom`: the chromosome
/// * `pos`: a genomic position
/// * `ref_seq`: the reference sequence observed at chrom:pos
/// * `variant_class`: "deletion", "insertion", or "substitution"
/// * `window_size`: number of base pairs to construct a centered window around the position
/// * `fasta_handle`: a file handle to the corresponding FASTA file for the sample's reference
///
/// returns: proportion GC content in [0., 1.] in a window centered around chrom:pos
fn calc_gc_proportion(
    chrom: &str,
    pos: u64,
    ref_seq: &str,
    variant_class: &str,
    window_size: u64,
    fasta_handle: &mut fasta::IndexedReader<File>,
) -> Result<f32, Box<dyn Error>> {
    // mutable position to adjust it based on the variant class
    let mut pos = pos;

    match variant_class {
        "deletion" => {
            // shift to the right to account for ref/alt alleles containing 1bp left context
            pos += 1;
        }
        "substitution" if ref_seq.len() > 2 => {
            // center window around middle of sequence (e.g. move 1 to the right for sequences of
            // length 3 or 4)
            pos += ((ref_seq.len() - 1) / 2) as u64;
        }
        _ => {}
    }

    // calculate the start and end of the symmetric window around the position
    let half_window = window_size / 2;
    let start = cmp::max(1, pos.saturating_sub(half_window));
    let end = pos + half_window;

    // fetch the sequence from the fasta file
    let mut sequence = Vec::new();
    fasta_handle.fetch(chrom, start - 1, end)?;
    fasta_handle.read(&mut sequence)?;

    // calculate GC proportion
    let gc_proportion = gc_content(sequence);

    Ok(gc_proportion)
}

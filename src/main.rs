fn main() -> Result<(), Box<dyn std::error::Error>> {
    use rust_htslib::{bam, bam::Read};
    use std::collections::HashSet;
    use std::convert::TryInto;
    use std::convert::TryFrom;
    let args: Vec<String> = std::env::args().collect();
    if args.len() - 1 != 5 {
        let mut usage = String::from("#arguments must be 5");
        usage.push_str("\nArg1: BAM file sorted by name");
        usage.push_str("\nArg2: Column header prefix for set1");
        usage.push_str("\nArg3: Set1 of reference names (comma separated)");
        usage.push_str("\nArg4: Column header prefix for set2");
        usage.push_str("\nArg5: Set2 of reference names (comma separated)");
        eprintln!("{}", usage);
        panic!("#arguments must be 5. observed: {}", args.len() - 1);
    }
    // Assuming the input bam is sorted by template name
    //let bam_path = &"mapped.alignmentset.merged.sorted_by_name.bam";
    let bam_path = &args[1];
    let set1_name = &args[2];
    let set1_refs: Vec<&str> = args[3].split(',').collect();
    let set2_name = &args[4];
    let set2_refs: Vec<&str> = args[5].split(',').collect();
    let mut bam = bam::Reader::from_path(bam_path)?;
    let header = bam.header();
    let num_tid: i32 = header.target_count().try_into().unwrap();
    eprintln!("INFO: target_count: {}", num_tid);
    (0 .. num_tid).for_each(|x| eprintln!("INFO: tid {} = {:?}", x, std::str::from_utf8(header.tid2name(x.try_into().unwrap())).unwrap()));
    let set1: HashSet<i32> = set1_refs.into_iter()
        .map(|ref_str| i32::try_from(header.tid(ref_str.as_bytes()).unwrap_or_else(|| panic!("Not found in the input bam: {}", ref_str))).unwrap()).collect();
    let set2: HashSet<i32> = set2_refs.into_iter()
        .map(|ref_str| i32::try_from(header.tid(ref_str.as_bytes()).unwrap_or_else(|| panic!("Not found in the input bam: {}", ref_str))).unwrap()).collect();
    set1.clone().into_iter().for_each(|x| eprintln!("INFO: set1: resolved tid is {} for {:?}", x, std::str::from_utf8(header.tid2name(x.try_into().unwrap())).unwrap()));
    set2.clone().into_iter().for_each(|x| eprintln!("INFO: set2: resolved tid is {} for {:?}", x, std::str::from_utf8(header.tid2name(x.try_into().unwrap())).unwrap()));
    let mut record = bam::Record::new();
    // Initial MapQ value for unobserved alignment
    // We may want to distinguish alignments with MapQ 0 from initial values
    let initial_mapq: u8 = 0;
    let mut set1_max = initial_mapq;
    let mut set2_max = initial_mapq;
    let mut set1_count: u32 = 0;
    let mut set2_count: u32 = 0;
    let mut num_valid: u64 = 0;
    let mut num_skipped_unmapped = 0;
    let mut num_skipped_outside_tid_sets = 0;
    let mut current_qname: String;
    let mut last_qname: String = "".to_string();
    let mut current_tid;
    let mut current_mapq: u8;
    use std::io::Write;
    let buffer_size: usize = 1024 * 1024; // 1M
    let stream = std::io::stdout();
    let mut buffer = std::io::BufWriter::with_capacity(buffer_size, stream.lock());
    // CSV header
    //println!("read_name,top_mapq_{},count_{},top_mapq_{},count_{},abs_diff_mapq", set1_name, set1_name, set2_name, set2_name);
    writeln!(buffer, "read_name,top_mapq_{},count_{},top_mapq_{},count_{},abs_diff_mapq", set1_name, set1_name, set2_name, set2_name)?;
    while bam.read(&mut record)? {
        current_qname = std::str::from_utf8(record.qname())?.to_string();
        current_tid = record.tid();
        current_mapq = record.mapq();
        if current_tid == -1 || current_mapq == 255 {
            num_skipped_unmapped += 1;
        } else if current_tid < num_tid {
            if current_qname != last_qname {
                if num_valid > 0 {
                    //println!("{},{},{},{},{},{}", last_qname, set1_max, set1_count, set2_max, set2_count, ((set1_max as i16) - (set2_max as i16)).abs());
                    writeln!(buffer, "{},{},{},{},{},{}", last_qname, set1_max, set1_count, set2_max, set2_count, ((set1_max as i16) - (set2_max as i16)).abs())?;
                }
                set1_max = initial_mapq;
                set2_max = initial_mapq;
                set1_count = 0;
                set2_count = 0;
            }
            if set1.contains(&current_tid) {
                set1_max = std::cmp::max(set1_max, current_mapq);
                set1_count += 1;
                num_valid += 1;
            } else if set2.contains(&current_tid) {
                set2_max = std::cmp::max(set2_max, current_mapq);
                set2_count += 1;
                num_valid += 1;
            } else {
                num_skipped_outside_tid_sets += 1;
            }
        } else {
            panic!("Unexpected tid: {}", current_tid);
        }
        last_qname = current_qname.clone();
        //println!("QNAME: {:?}, TID: {:?}, MAPQ: {:?}", current_qname, current_tid, current_mapq);
    }
    // Write the result for the last block
    if num_valid > 0 {
        //println!("{},{},{},{},{},{}", last_qname, set1_max, set1_count, set2_max, set2_count, ((set1_max as i16) - (set2_max as i16)).abs());
        writeln!(buffer, "{},{},{},{},{},{}", last_qname, set1_max, set1_count, set2_max, set2_count, ((set1_max as i16) - (set2_max as i16)).abs())?;
    }
    eprintln!("INFO: # valid alignments: {}", num_valid);
    eprintln!("INFO: # unmapped reads: {}", num_skipped_unmapped);
    eprintln!("INFO: # alignments outside specified sets: {}", num_skipped_outside_tid_sets);
    buffer.flush()?;
    Ok(())
}

use fixed::types::U0F64;
use num::{Bounded, FromPrimitive, ToPrimitive};
use std::fs::File;
use std::io::Write;
use std::{env, process};

fn build_cordic_atan_table(num_entries: u8) -> Vec<u64> {
    let mut table = Vec::with_capacity(num_entries as usize);
    let mut angle = 1.0f64;

    for _ in 0..num_entries {
        let atan = angle.atan();
        let ival = U0F64::from_num(atan);
        table.push(ival.to_bits());
        angle = angle * 0.5;
    }

    table
}

// NOTE: we generate exp(x) - 1.0 instead of exp(x) directly
// so that the results can fit into a U0F64 like all our other
// CORDIC lookup tables.
fn build_cordic_exp_minus_one_table(num_entries: u8) -> Vec<u64> {
    let mut table = Vec::with_capacity(num_entries as usize);
    let mut angle = 0.5f64;

    for _ in 0..num_entries {
        let exp = angle.exp();
        let ival = U0F64::from_num(exp - 1.0);
        table.push(ival.to_bits());
        angle = angle * 0.5;
    }

    table
}

fn build_cordic_cumprod_table(num_entries: u8) -> Vec<u64> {
    let mut table = Vec::with_capacity(num_entries as usize);
    let mut cumprod = 1.0f64;
    let mut pow = 1.0f64;

    for _ in 0..num_entries {
        cumprod /= (1.0 + pow).sqrt();
        let ival = U0F64::from_num(cumprod);
        table.push(ival.to_bits());
        pow *= 0.25;
    }

    table
}

fn build_sin_table<T: Bounded + ToPrimitive + FromPrimitive + std::fmt::Display>(
    num_bits: u8,
) -> Vec<T> {
    let num: u32 = 1 << num_bits;
    let mut table = Vec::with_capacity((num + 1) as usize);
    for i in 0..(num + 1) {
        let angle = (std::f64::consts::FRAC_PI_2 * f64::from(i)) / f64::from(num);
        let val = (angle.sin() * T::max_value().to_f64().unwrap()).round();
        let ival = T::from_f64(val).unwrap();
        table.push(ival);
    }
    table
}

fn build_asin_table<T: Bounded + ToPrimitive + FromPrimitive + std::fmt::Display>(
    num_bits: u8,
) -> Vec<T> {
    let num: u32 = 1 << num_bits;
    let mut table = Vec::with_capacity((num + 1) as usize);
    for i in 0..num {
        let val = f64::from(i) / f64::from(num);
        let angle =
            ((val.asin() * T::max_value().to_f64().unwrap()) / std::f64::consts::FRAC_PI_2).round();
        let ival = T::from_f64(angle).unwrap();
        table.push(ival);
    }
    table.push(T::max_value());
    table
}

fn vec_u16_to_le_bytes(data: Vec<u16>) -> Vec<u8> {
    let mut bytes = vec![];
    data.into_iter()
        .for_each(|d| bytes.extend(&d.to_le_bytes()));
    bytes
}

fn vec_u64_to_le_bytes(data: Vec<u64>) -> Vec<u8> {
    let mut bytes = vec![];
    data.into_iter()
        .for_each(|d| bytes.extend(&d.to_le_bytes()));
    bytes
}

fn write_table(table_filepath: &str, num_bits: Option<u8>, data: Vec<u8>) {
    let mut table_file = File::create(table_filepath).expect("unable to open file");
    if let Some(num_bits) = num_bits {
        table_file
            .write_all(&[num_bits])
            .expect("unable to write to file");
    }
    table_file
        .write_all(&data)
        .expect("unable to write to file");
}

fn print_usage(error: Option<&str>) {
    if let Some(error) = error {
        eprintln!("Error: {}", error);
    }
    println!(
        "usage: table-builder {{sin|asin|cordic-atan|cordic-exp-minus-one|cordic-cumprod}} num_bits table_filepath"
    );
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() == 1 {
        print_usage(None);
        process::exit(0);
    }
    if args.len() != 4 {
        print_usage(Some("Incorrect number of parameters"));
        process::exit(1);
    }
    let num = args[2].parse::<u8>().expect("unable to parse num_bits");
    match args[1].as_str() {
        "sin" => {
            let table = build_sin_table::<u16>(num as u8);
            let data = vec_u16_to_le_bytes(table);
            write_table(args[3].as_str(), Some(num), data);
        }
        "asin" => {
            let table = build_asin_table::<u16>(num as u8);
            let data = vec_u16_to_le_bytes(table);
            write_table(args[3].as_str(), Some(num), data);
        }
        "cordic-atan" => {
            let table = build_cordic_atan_table(num);
            let data = vec_u64_to_le_bytes(table);
            write_table(args[3].as_str(), None, data);
        }
        "cordic-exp-minus-one" => {
            let table = build_cordic_exp_minus_one_table(num);
            let data = vec_u64_to_le_bytes(table);
            write_table(args[3].as_str(), None, data);
        }
        "cordic-cumprod" => {
            let table = build_cordic_cumprod_table(num);
            let data = vec_u64_to_le_bytes(table);
            write_table(args[3].as_str(), None, data);
        }
        _ => {
            print_usage(Some("Invalid table type"));
            process::exit(1)
        }
    }
}

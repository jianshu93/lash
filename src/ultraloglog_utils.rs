use std::io::{BufRead, BufReader, Write, BufWriter};
use std::fs::File;
use std::path::Path;
use std::str::FromStr;
use std::fmt;

pub fn load<T>(filepath: &str) -> Vec<T>
where
    T: FromStr + fmt::Debug,
{
    let reader = BufReader::new(File::open(filepath).unwrap());
    let mut nums = Vec::with_capacity(10000);

    for line in reader.lines() {
        nums.push(
            line.unwrap()
                .parse::<T>()
                .map_err(|_| "Parsing line failed")
                .unwrap(),
        );
    }

    nums
}

// Saves values to a file
pub fn save<T>(values: &Vec<T>, filename: &str, output: &str)
where
    T: fmt::Display,
{
    let mut writer = BufWriter::new(File::create(Path::new(output).join(filename)).unwrap());

    for val in values {
        write!(writer, "{}\n", val).unwrap();
    }

    writer.flush().unwrap();
}

pub struct Estimation(pub u64, pub u64, pub String);

impl fmt::Display for Estimation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} {} {}", self.0, self.1, self.2)
    }
}
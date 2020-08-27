extern crate byteorder;
extern crate cpython;
extern crate flate2;

use std::error::Error;
use std::fs::File;
use std::io;
use std::io::{BufReader, Read};
use std::mem;
use std::path::Path;

use byteorder::{ByteOrder, LittleEndian};
use flate2::read::GzDecoder;

// exposing the library as a cython for vpolo EDS reading
use cpython::{py_fn, py_module_initializer, PyResult, Python};

py_module_initializer!(sce, |py, m| {
    m.add(py, "__doc__", "This module is implemented in Rust.")?;
    m.add(
        py,
        "read_quants",
        py_fn!(py, read_eds(filename: &str, num_rows: u64, num_cols: u64)),
    )?;
    Ok(())
});

fn read_eds(_py: Python, val: &str, num_rows: u64, num_cols: u64) -> PyResult<(Vec<usize>, Vec<usize>, Vec<MatValT>)> {
    let fpath = Path::new(val);
    let mat = reader(fpath, num_rows as usize, num_cols as usize).expect("can't read the EDS file");

    Ok((mat.i, mat.j, mat.k))
}

pub type MatValT = f32;
pub struct SpMatrix {
    pub i: Vec<usize>,
    pub j: Vec<usize>,
    pub k: Vec<MatValT>,
}

fn get_reserved_spaces(
    num_bit_vecs: usize,
    num_rows: usize,
    mut file: GzDecoder<BufReader<File>>,
) -> Result<Vec<usize>, Box<dyn Error>> {
    let mut bit_vec_lengths: Vec<usize> = vec![0; num_rows + 1];
    let mut bit_vec = vec![0; num_bit_vecs];
    let mut running_sum = 0;

    for i in 0..num_rows {
        file.read_exact(&mut bit_vec[..])?;
        let mut num_ones = 0;
        for bits in bit_vec.iter() {
            num_ones += bits.count_ones() as usize;
        }

        running_sum += num_ones;
        bit_vec_lengths[i + 1] = running_sum;

        // no seek command yet
        // copied from https://github.com/rust-lang/rust/issues/53294#issue-349837288
        io::copy(
            &mut file
                .by_ref()
                .take((num_ones * mem::size_of::<MatValT>()) as u64),
            &mut io::sink(),
        )?;
        //file.seek(SeekFrom::Current((num_ones * 4) as i64))?;
    }

    Ok(bit_vec_lengths)
}

// reads the EDS format single cell matrix from the given path
pub fn reader(
    file_path: &Path,
    num_rows: usize,
    num_cols: usize,
) -> Result<SpMatrix, Box<dyn Error>> {
    // reading the matrix
    let file_handle = File::open(file_path.to_owned())?;
    let buffered = BufReader::new(file_handle);
    let file = GzDecoder::new(buffered);

    let num_bit_vecs: usize = (num_cols + 7) / 8;
    let bit_vector_lengths = get_reserved_spaces(num_bit_vecs, num_rows, file)?;

    let total_nnz = bit_vector_lengths[num_rows];
    let mut data: Vec<MatValT> = vec![0.0; total_nnz];
    let mut indices: Vec<usize> = vec![0; total_nnz];

    let file_handle = File::open(file_path)?;
    let buffered = BufReader::new(file_handle);
    let mut file = GzDecoder::new(buffered);

    let mut global_pointer = 0;
    let mut bit_vec = vec![0; num_bit_vecs];
    for i in 0..num_rows {
        file.read_exact(&mut bit_vec[..])?;
        let num_ones = bit_vector_lengths[i + 1] - bit_vector_lengths[i];

        let mut one_validator = 0;
        for (j, flag) in bit_vec.iter().enumerate() {
            if *flag != 0 {
                for (i, bit_id) in format!("{:8b}", flag).chars().enumerate() {
                    if let '1' = bit_id {
                        let offset = i + (8 * j);
                        indices[global_pointer + one_validator] = offset;

                        one_validator += 1;
                    };
                }
            }
        }
        assert_eq!(num_ones, one_validator);

        let mut expression: Vec<u8> = vec![0; mem::size_of::<MatValT>() * (num_ones as usize)];
        let mut float_buffer: Vec<MatValT> = vec![0.0; num_ones as usize];
        file.read_exact(&mut expression[..])?;

        // NOTE: if we change MatValT, double check below line
        LittleEndian::read_f32_into(&expression, &mut float_buffer);

        for value in float_buffer {
            data[global_pointer] = value as MatValT;
            global_pointer += 1;
        }
    }

    assert_eq!(global_pointer, total_nnz);
    let matrix = SpMatrix{i: bit_vector_lengths, j: indices, k: data};
    Ok(matrix)
}

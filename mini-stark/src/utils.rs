use crate::merkle::MerkleTree;
use crate::Column;
use ark_poly::domain::Radix2EvaluationDomain;
use ark_poly::EvaluationDomain;
use digest::Digest;
use digest::Output;
use fast_poly::allocator::PageAlignedAllocator;
use fast_poly::plan::Fft;
use fast_poly::plan::Ifft;
use fast_poly::GpuField;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::ops::Index;
use std::ops::IndexMut;
use std::time::Instant;

/// Matrix is an array of columns.
pub struct Matrix<F>(Vec<Vec<F, PageAlignedAllocator>>);

impl<F: GpuField> Matrix<F> {
    pub fn new(cols: Vec<Vec<F, PageAlignedAllocator>>) -> Self {
        Matrix(cols)
    }

    // TODO: perhaps bring naming of rows and cols in line with
    // how the trace is names i.e. len and width.
    pub fn num_rows(&self) -> usize {
        if self.0.is_empty() {
            return 0;
        }
        // Check all columns have the same length
        let expected_len = self.0[0].len();
        assert!(self.0.iter().all(|col| col.len() == expected_len));
        expected_len
    }

    pub fn num_cols(&self) -> usize {
        self.0.len()
    }

    pub fn is_empty(&self) -> bool {
        self.num_rows() == 0
    }

    pub fn interpolate_columns(&self) -> Self {
        let domain = Radix2EvaluationDomain::<F>::new(self.num_rows()).unwrap();
        let mut ifft = Ifft::from(domain);
        let columns = self
            .0
            .iter()
            .map(|evaluations| {
                let mut poly = evaluations.to_vec_in(PageAlignedAllocator);
                ifft.encode(&mut poly);
                poly
            })
            .collect();

        let start = Instant::now();
        ifft.execute();
        println!("THE Completion: {:?}", start.elapsed());

        Matrix::new(columns)
    }

    pub fn evaluate(&self, domain: Radix2EvaluationDomain<F>) -> Self {
        let mut fft = {
            let _timer = Timer::new("FFT eval creation");
            Fft::from(domain)
        };
        let columns = {
            let _timer = Timer::new("Actual evaluations");
            self.0
                .iter()
                .map(|poly| {
                    let mut evaluations = poly.to_vec_in(PageAlignedAllocator);
                    fft.encode(&mut evaluations);
                    evaluations
                })
                .collect()
        };

        let start = Instant::now();
        fft.execute();
        println!("THE Completion: {:?}", start.elapsed());

        Matrix::new(columns)
    }

    pub fn commit_to_rows<D: Digest>(&self) -> MerkleTree<D> {
        let num_rows = self.num_rows();
        let num_cols = self.num_cols();
        let mut row_hashes = Vec::with_capacity(num_rows);
        for row in 0..num_rows {
            // TODO: move this outside the for loop and use .reset()
            let mut hasher = D::new();
            for col in 0..num_cols {
                let mut bytes = Vec::new();
                self.0[col][row].serialize_compressed(&mut bytes);
                hasher.update(&bytes);
            }
            row_hashes.push(hasher.finalize())
        }
        MerkleTree::new(row_hashes).expect("failed to construct Merkle tree")
    }
}

impl<F: GpuField, C: Column> Index<C> for Matrix<F> {
    type Output = Vec<F, PageAlignedAllocator>;

    fn index(&self, col: C) -> &Self::Output {
        &self.0[col.index()]
    }
}

impl<F: GpuField, C: Column> IndexMut<C> for Matrix<F> {
    fn index_mut(&mut self, col: C) -> &mut Self::Output {
        &mut self.0[col.index()]
    }
}

pub struct Timer<'a> {
    name: &'a str,
    start: Instant,
}

impl<'a> Timer<'a> {
    pub fn new(name: &'a str) -> Timer<'a> {
        let start = Instant::now();
        Timer { name, start }
    }
}

impl<'a> Drop for Timer<'a> {
    fn drop(&mut self) {
        println!("{} in {:?}", self.name, self.start.elapsed());
    }
}

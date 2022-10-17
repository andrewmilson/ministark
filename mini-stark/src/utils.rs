use crate::merkle::MerkleTree;
use crate::Column;
use ark_ff::Zero;
use ark_poly::domain::Radix2EvaluationDomain;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_poly::Polynomial;
use digest::Digest;
use fast_poly::allocator::PageAlignedAllocator;
use fast_poly::plan::GpuFft;
use fast_poly::plan::GpuIfft;
use fast_poly::plan::PLANNER;
use fast_poly::stage::AddAssignStage;
use fast_poly::utils::buffer_mut_no_copy;
use fast_poly::utils::buffer_no_copy;
use fast_poly::GpuField;
use fast_poly::GpuVec;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::ops::Deref;
use std::ops::DerefMut;
use std::ops::Index;
use std::ops::IndexMut;
use std::time::Instant;

/// Matrix is an array of columns.
pub struct Matrix<F>(pub Vec<GpuVec<F>>);

impl<F: GpuField> Matrix<F> {
    pub fn new(cols: Vec<GpuVec<F>>) -> Self {
        Matrix(cols)
    }

    pub fn from_rows(rows: Vec<Vec<F>>) -> Self {
        let num_rows = rows.len();
        let num_cols = rows.first().map(|first| first.len()).unwrap_or(0);
        let mut cols = (0..num_cols)
            .map(|_| Vec::with_capacity_in(num_rows, PageAlignedAllocator))
            .collect::<Vec<GpuVec<F>>>();
        // TODO: parallelise
        for row in rows {
            debug_assert_eq!(row.len(), num_cols);
            for (col, value) in cols.iter_mut().zip(row) {
                col.push(value)
            }
        }
        Matrix::new(cols)
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

    pub fn append(&mut self, other: Matrix<F>) {
        for col in other.0 {
            self.0.push(col)
        }
    }

    pub fn join(mut matrices: Vec<Matrix<F>>) -> Matrix<F> {
        let mut accumulator = Vec::new();
        for matrix in &mut matrices {
            accumulator.append(matrix)
        }
        Matrix::new(accumulator)
    }

    pub fn num_cols(&self) -> usize {
        self.0.len()
    }

    pub fn is_empty(&self) -> bool {
        self.num_rows() == 0
    }

    pub fn interpolate_columns(&self, domain: Radix2EvaluationDomain<F>) -> Self {
        let mut ifft = GpuIfft::from(domain);
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
            GpuFft::from(domain)
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

    /// Sums columns into a single column matrix
    pub fn sum_columns(&self) -> Matrix<F> {
        let n = self.num_rows();
        let mut accumulator = Vec::with_capacity_in(n, PageAlignedAllocator);
        accumulator.resize(n, F::zero());

        if !self.num_cols().is_zero() {
            // TODO: could improve
            let library = &PLANNER.library;
            let command_queue = &PLANNER.command_queue;
            let device = command_queue.device();
            let command_buffer = command_queue.new_command_buffer();
            let mut accumulator_buffer = buffer_mut_no_copy(device, &mut accumulator);
            let adder = AddAssignStage::<F>::new(library, n);
            for column in self.0.iter() {
                let column_buffer = buffer_no_copy(command_queue.device(), column);
                adder.encode(command_buffer, &mut accumulator_buffer, &column_buffer);
            }
            command_buffer.commit();
            command_buffer.wait_until_completed();
        }

        Matrix::new(vec![accumulator])
    }

    pub fn commit_to_rows<D: Digest>(&self) -> MerkleTree<D> {
        let num_rows = self.num_rows();
        let num_cols = self.num_cols();

        let mut row_hashes = Vec::with_capacity(num_rows);
        row_hashes.resize(num_rows, Default::default());
        ark_std::cfg_iter_mut!(row_hashes)
            .enumerate()
            .for_each(|(row, row_hash)| {
                let mut hasher = D::new();
                // TODO: bad for single column matricies
                for col in 0..num_cols {
                    let mut bytes = Vec::new();
                    self.0[col][row].serialize_compressed(&mut bytes).unwrap();
                    hasher.update(&bytes);
                }
                *row_hash = hasher.finalize();
            });
        MerkleTree::new(row_hashes).expect("failed to construct Merkle tree")
    }

    pub fn evaluate_at(&self, x: F) -> Vec<F> {
        let mut evaluations = Vec::new();
        for col in &self.0 {
            // TODO: perf
            let polynomial = DensePolynomial::from_coefficients_slice(col);
            evaluations.push(polynomial.evaluate(&x));
        }
        evaluations
    }

    pub fn get_row(&self, row: usize) -> Option<Vec<F>> {
        if row < self.num_rows() {
            Some(self.iter().map(|col| col[row]).collect())
        } else {
            None
        }
    }

    pub fn rows(&self) -> Vec<Vec<F>> {
        (0..self.num_rows())
            .map(|row| self.get_row(row).unwrap())
            .collect()
    }

    pub fn column_degrees(&self) -> Vec<usize> {
        self.0
            .iter()
            .map(|col| {
                let mut num_leading_zeros = 0;
                while col.last().map_or(false, |c| c.is_zero()) {
                    num_leading_zeros += 1;
                }
                col.len() - num_leading_zeros - 1
            })
            .collect()
    }
}

impl<F: GpuField> Clone for Matrix<F> {
    fn clone(&self) -> Self {
        Self(
            self.0
                .iter()
                .map(|col| col.to_vec_in(PageAlignedAllocator))
                .collect(),
        )
    }
}

impl<F: GpuField> DerefMut for Matrix<F> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl<F: GpuField> Deref for Matrix<F> {
    type Target = Vec<GpuVec<F>>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<F: GpuField, C: Column> Index<C> for Matrix<F> {
    type Output = GpuVec<F>;

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

pub fn interleave<T: Copy + Send + Sync + Default, const N: usize>(source: &[T]) -> Vec<[T; N]> {
    let n = source.len() / N;
    let mut res = vec![[T::default(); N]; n];
    ark_std::cfg_iter_mut!(res)
        .enumerate()
        .for_each(|(i, element)| {
            for j in 0..N {
                element[j] = source[i + j * n]
            }
        });
    res
}

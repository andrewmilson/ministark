use crate::Column;
use fast_poly::allocator::PageAlignedAllocator;
use fast_poly::GpuField;
use std::ops::Index;
use std::ops::IndexMut;

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

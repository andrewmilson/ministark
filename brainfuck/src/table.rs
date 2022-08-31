use std::marker::PhantomData;

pub struct Table<E> {
    base_width: usize,
    full_width: usize,
    length: usize,
    num_randomizers: usize,
    height: usize,
    _phantom: PhantomData<E>,
}

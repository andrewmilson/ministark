// Implementation is adapted from RationalExpression in https://github.com/0xProject/OpenZKP
#![allow(clippy::arc_with_non_send_sync)]

use alloc::collections::BTreeMap;
use ark_ff::One;
use ark_std::Zero;
use core::cmp::Ordering;
use core::hash::Hash;
use core::iter::Product;
use core::iter::Sum;
use core::ops::Add;
use core::ops::AddAssign;
use core::ops::Div;
use core::ops::DivAssign;
use core::ops::Mul;
use core::ops::MulAssign;
use core::ops::Neg;
use core::ops::Sub;
use core::ops::SubAssign;
use num_traits::Pow;
use std::cell::RefCell;
use std::collections::hash_map::DefaultHasher;
use std::hash::Hasher;
use std::ptr::addr_of;
use std::rc::Rc;
use std::sync::Arc;
use std::sync::RwLock;

pub type P<T> = Arc<RwLock<T>>;

/// Expression
#[derive(Clone)]
pub enum Expr<T> {
    Leaf(T),
    Neg(P<Expr<T>>),
    Add(P<Expr<T>>, P<Expr<T>>),
    Mul(P<Expr<T>>, P<Expr<T>>),
    Div(P<Expr<T>>, P<Expr<T>>),
    Pow(P<Expr<T>>, usize),
}

// unsafe impl<T: Send> Send for Expr<T> {}
// unsafe impl<T: Sync> Sync for Expr<T> {}

impl<T> Expr<T> {
    // Adapted from https://github.com/0xProject/OpenZKP
    /// Applies a bottom-up traversal.
    pub fn traverse(&self, f: &mut impl FnMut(&Self)) {
        use Expr::*;
        match self {
            // Tree types are recursed first
            Neg(a) | Pow(a, _) => a.as_ref().read().unwrap().traverse(f),
            Add(a, b) | Mul(a, b) | Div(a, b) => {
                a.as_ref().read().unwrap().traverse(f);
                b.as_ref().read().unwrap().traverse(f);
            }
            Leaf(_) => {}
        }

        f(self);
    }

    // Adapted from https://github.com/0xProject/OpenZKP
    /// Applies a mapped bottom-up traversal.
    /// The function applies to each node after application to its descendants
    pub fn map(&self, f: &mut impl FnMut(Self) -> Self) -> Self
    where
        T: Clone,
    {
        use Expr::*;
        // TODO: consider function argument for leaf nodes
        // TODO: why can't the copiler do this as a param
        let res = match self {
            // Tree types are recursed first
            Add(a, b) => Add(
                Arc::new(RwLock::new(a.as_ref().read().unwrap().map(f))),
                Arc::new(RwLock::new(b.as_ref().read().unwrap().map(f))),
            ),
            Neg(a) => Neg(Arc::new(RwLock::new(a.as_ref().read().unwrap().map(f)))),
            Mul(a, b) => Mul(
                Arc::new(RwLock::new(a.as_ref().read().unwrap().map(f))),
                Arc::new(RwLock::new(b.as_ref().read().unwrap().map(f))),
            ),
            Div(a, b) => Div(
                Arc::new(RwLock::new(a.as_ref().read().unwrap().map(f))),
                Arc::new(RwLock::new(b.as_ref().read().unwrap().map(f))),
            ),
            Pow(a, e) => Pow(Arc::new(RwLock::new(a.as_ref().read().unwrap().map(f))), *e),

            // Leaf types are mapped as is.
            Leaf(v) => Leaf(v.clone()),
        };

        f(res)
    }

    /// Applies a bottom-up traversal.
    /// The closure is given mutable access to the nodes.
    pub fn traverse_mut(&mut self, f: &mut impl FnMut(&mut Self)) {
        use Expr::*;
        match self {
            // Tree types are recursed first
            Add(a, b) | Mul(a, b) | Div(a, b) => {
                a.write().unwrap().traverse_mut(f);
                b.write().unwrap().traverse_mut(f);
            }
            Neg(a) | Pow(a, _) => a.write().unwrap().traverse_mut(f),
            Leaf(_) => {}
        }

        f(self);
    }

    // TODO: change this to self
    fn _map_leaves<U>(
        this: &P<Self>,
        sl: &mut BTreeMap<*const T, P<Expr<U>>>,
        sn: &mut BTreeMap<*const Self, P<Expr<U>>>,
        f: &mut impl FnMut(&T) -> U,
    ) -> P<Expr<U>>
    where
        T: Ord,
    {
        // The lint's suggestion "Using Option::map_or_else" doesn't work because
        // `.insert` requires mutable access to seen (which is already borrowed).
        #[allow(clippy::option_if_let_else)]
        if let Some(node) = sn.get(&addr_of!(*this.read().unwrap())) {
            Arc::clone(node)
        } else {
            use Expr::*;
            let res = Arc::new(RwLock::new(match &*this.read().unwrap() {
                Leaf(l) => {
                    if let Some(leaf) = sl.get(&addr_of!(*l)) {
                        return Arc::clone(leaf);
                    }

                    let res = Arc::new(RwLock::new(Leaf(f(l))));
                    sl.insert(addr_of!(*l), Arc::clone(&res));
                    sn.insert(addr_of!(*this.read().unwrap()), Arc::clone(&res));
                    return res;
                }
                Add(a, b) => Add(
                    Self::_map_leaves(a, sl, sn, f),
                    Self::_map_leaves(b, sl, sn, f),
                ),
                Neg(a) => Neg(Self::_map_leaves(a, sl, sn, f)),
                Pow(a, e) => Pow(Self::_map_leaves(a, sl, sn, f), *e),
                Mul(a, b) => Mul(
                    Self::_map_leaves(a, sl, sn, f),
                    Self::_map_leaves(b, sl, sn, f),
                ),
                Div(a, b) => Div(
                    Self::_map_leaves(a, sl, sn, f),
                    Self::_map_leaves(b, sl, sn, f),
                ),
            }));
            sn.insert(addr_of!(*this.read().unwrap()), Arc::clone(&res));
            res
        }
    }

    // Maps leaves and retains internal structure
    #[allow(clippy::many_single_char_names)]
    pub fn map_leaves<U>(&self, f: &mut impl FnMut(&T) -> U) -> Expr<U>
    where
        T: Ord,
    {
        use Expr::*;
        let mut leafs_map = BTreeMap::new();
        let mut nodes_map = BTreeMap::new();
        let l = &mut leafs_map;
        let m = &mut nodes_map;
        match self {
            Leaf(v) => Leaf(f(v)),
            Add(a, b) => Add(Self::_map_leaves(a, l, m, f), Self::_map_leaves(b, l, m, f)),
            Neg(a) => Neg(Self::_map_leaves(a, l, m, f)),
            Mul(a, b) => Mul(Self::_map_leaves(a, l, m, f), Self::_map_leaves(b, l, m, f)),
            Div(a, b) => Div(Self::_map_leaves(a, l, m, f), Self::_map_leaves(b, l, m, f)),
            Pow(a, e) => Pow(Self::_map_leaves(a, l, m, f), *e),
        }
    }

    /// Inspired by <https://neptune.cash/learn/speed-up-stark-provers-with-multicircuits/>
    /// Runtime: O(n log n) where n is the number of edges TODO: check
    #[allow(clippy::too_many_lines)]
    pub fn reuse_shared_nodes(&self) -> Self
    where
        T: Ord + Clone + Hash,
    {
        type Id = u64;
        type SeenSet<T> = Rc<RefCell<BTreeMap<Id, P<Expr<T>>>>>;

        struct IdNode<T> {
            node: P<Expr<T>>,
            seen: SeenSet<T>,
            id: Id,
        }

        impl<T: Hash + Clone> IdNode<T> {
            fn new_leaf(leaf: &T, seen: SeenSet<T>) -> Self {
                // `id` is the hash of the leaf
                let mut hasher = DefaultHasher::new();
                ("leaf", leaf).hash(&mut hasher);
                let id = hasher.finish();

                // can't use if_let_else because it causes borrow errors
                #[allow(clippy::option_if_let_else)]
                let node = if seen.borrow().contains_key(&id) {
                    // we have encountered this leaf
                    Arc::clone(seen.borrow().get(&id).unwrap())
                } else {
                    // we haven't encountered this leaf, create a new entry
                    let node = Arc::new(RwLock::new(Expr::Leaf(leaf.clone())));
                    seen.borrow_mut().insert(id, Arc::clone(&node));
                    node
                };

                Self { node, seen, id }
            }
        }

        impl<T> Add for IdNode<T> {
            type Output = Self;

            fn add(self, rhs: Self) -> Self::Output {
                let seen = self.seen;

                let mut hasher = DefaultHasher::new();
                ("add", self.id, rhs.id).hash(&mut hasher);
                let id = hasher.finish();

                // can't use if_let_else because it causes borrow errors
                #[allow(clippy::option_if_let_else)]
                let node = if seen.borrow().contains_key(&id) {
                    // we have encountered this node
                    Arc::clone(seen.borrow().get(&id).unwrap())
                } else {
                    // we haven't encountered this node, create a new entry
                    let node = Arc::new(RwLock::new(Expr::Add(self.node, rhs.node)));
                    seen.borrow_mut().insert(id, Arc::clone(&node));
                    node
                };

                Self { node, seen, id }
            }
        }

        impl<T> Mul for IdNode<T> {
            type Output = Self;

            fn mul(self, rhs: Self) -> Self::Output {
                let seen = self.seen;

                let mut hasher = DefaultHasher::new();
                ("mul", self.id, rhs.id).hash(&mut hasher);
                let id = hasher.finish();

                // can't use if_let_else because it causes borrow errors
                #[allow(clippy::option_if_let_else)]
                let node = if seen.borrow().contains_key(&id) {
                    // we have encountered this leaf
                    Arc::clone(seen.borrow().get(&id).unwrap())
                } else {
                    // we haven't encountered this node, create a new entry
                    let node = Arc::new(RwLock::new(Expr::Mul(self.node, rhs.node)));
                    seen.borrow_mut().insert(id, Arc::clone(&node));
                    node
                };

                Self { node, seen, id }
            }
        }

        impl<T> Div for IdNode<T> {
            type Output = Self;

            fn div(self, rhs: Self) -> Self::Output {
                let seen = self.seen;

                let mut hasher = DefaultHasher::new();
                ("div", self.id, rhs.id).hash(&mut hasher);
                let id = hasher.finish();

                // can't use if_let_else because it causes borrow errors
                #[allow(clippy::option_if_let_else)]
                let node = if seen.borrow().contains_key(&id) {
                    // we have encountered this leaf
                    Arc::clone(seen.borrow().get(&id).unwrap())
                } else {
                    // we haven't encountered this node, create a new entry
                    let node = Arc::new(RwLock::new(Expr::Div(self.node, rhs.node)));
                    seen.borrow_mut().insert(id, Arc::clone(&node));
                    node
                };

                Self { node, seen, id }
            }
        }

        impl<T> Neg for IdNode<T> {
            type Output = Self;

            fn neg(self) -> Self::Output {
                let seen = self.seen;

                let mut hasher = DefaultHasher::new();
                ("neg", self.id).hash(&mut hasher);
                let id = hasher.finish();

                // can't use if_let_else because it causes borrow errors
                #[allow(clippy::option_if_let_else)]
                let node = if seen.borrow().contains_key(&id) {
                    // we have encountered this leaf
                    Arc::clone(seen.borrow().get(&id).unwrap())
                } else {
                    // we haven't encountered this node, create a new entry
                    let node = Arc::new(RwLock::new(Expr::Neg(self.node)));
                    seen.borrow_mut().insert(id, Arc::clone(&node));
                    node
                };

                Self { node, seen, id }
            }
        }

        impl<T> Pow<usize> for IdNode<T> {
            type Output = Self;

            fn pow(self, exp: usize) -> Self::Output {
                let seen = self.seen;

                let mut hasher = DefaultHasher::new();
                ("pow", self.id, exp).hash(&mut hasher);
                let id = hasher.finish();

                // can't use if_let_else because it causes borrow errors
                #[allow(clippy::option_if_let_else)]
                let node = if seen.borrow().contains_key(&id) {
                    // we have encountered this leaf
                    Arc::clone(seen.borrow().get(&id).unwrap())
                } else {
                    // we haven't encountered this node, create a new entry
                    let node = Arc::new(RwLock::new(Expr::Pow(self.node, exp)));
                    seen.borrow_mut().insert(id, Arc::clone(&node));
                    node
                };

                Self { node, seen, id }
            }
        }

        let seen = Rc::new(RefCell::new(BTreeMap::new()));
        let res = self.eval(&mut |leaf| IdNode::new_leaf(leaf, Rc::clone(&seen)));
        // Drop references
        drop((seen, res.seen));
        Arc::into_inner(res.node).unwrap().into_inner().unwrap()
    }

    // Adapted from https://github.com/0xProject/OpenZKP
    // NOTE: evaluates the expression as a tree not a DAG
    /// Evaluates an expression bottom up as a tree
    /// Use `eval_graph` to evaluate and expression as a graph
    pub fn eval<U>(&self, f: &mut impl FnMut(&T) -> U) -> U
    where
        U: Add<Output = U>
            + Neg<Output = U>
            + Div<Output = U>
            + Mul<Output = U>
            + Pow<usize, Output = U>,
    {
        use Expr::*;
        // #[cfg(feature = "std")]
        // println!("brr");
        match self {
            Leaf(a) => f(a),
            Add(a, b) => a.read().unwrap().eval(f) + b.read().unwrap().eval(f),
            Neg(a) => -a.read().unwrap().eval(f),
            Mul(a, b) => a.read().unwrap().eval(f) * b.read().unwrap().eval(f),
            Div(a, b) => a.read().unwrap().eval(f) / b.read().unwrap().eval(f),
            Pow(a, e) => a.read().unwrap().eval(f).pow(*e),
        }
    }

    fn _graph_eval<U>(this: &P<Expr<Option<U>>>) -> U
    where
        U: Clone
            + Add<Output = U>
            + Neg<Output = U>
            + Div<Output = U>
            + Mul<Output = U>
            + Pow<usize, Output = U>,
    {
        use Expr::*;
        let mut node = this.write().unwrap();
        let res = match core::mem::take(&mut *node) {
            Leaf(v) => v.unwrap(),
            Neg(a) => -Self::_graph_eval(&a),
            Add(a, b) => Self::_graph_eval(&a) + Self::_graph_eval(&b),
            Mul(a, b) => Self::_graph_eval(&a) * Self::_graph_eval(&b),
            Div(a, b) => Self::_graph_eval(&a) / Self::_graph_eval(&b),
            Pow(a, e) => Self::_graph_eval(&a).pow(e),
        };
        *node = Leaf(Some(res.clone()));
        res
    }

    /// Evaluates an expression graph bottom up
    /// Intermediate results are cached to prevent re-evaluation
    pub fn graph_eval<U>(&self, f: &mut impl FnMut(&T) -> U) -> U
    where
        T: Ord + Copy,
        U: Clone
            + Add<Output = U>
            + Neg<Output = U>
            + Div<Output = U>
            + Mul<Output = U>
            + Pow<usize, Output = U>,
    {
        use Expr::*;
        let res = self.map_leaves(&mut |l| Some(f(l)));
        match res {
            Leaf(a) => a.unwrap(),
            Add(a, b) => Self::_graph_eval(&a) + Self::_graph_eval(&b),
            Neg(a) => -Self::_graph_eval(&a),
            Mul(a, b) => Self::_graph_eval(&a) * Self::_graph_eval(&b),
            Div(a, b) => Self::_graph_eval(&a) / Self::_graph_eval(&b),
            Pow(a, e) => Self::_graph_eval(&a).pow(e),
        }
    }
}

impl<T: Default> Default for Expr<T> {
    fn default() -> Self {
        Self::Leaf(T::default())
    }
}

impl<T: Eq> Eq for Expr<T> {}

impl<T: PartialEq> PartialEq for Expr<T> {
    fn eq(&self, other: &Self) -> bool {
        use Expr::*;
        match (self, other) {
            (Add(l0, l1), Add(r0, r1))
            | (Mul(l0, l1), Mul(r0, r1))
            | (Div(l0, l1), Div(r0, r1)) => {
                *l0.read().unwrap() == *r0.read().unwrap()
                    && *l1.read().unwrap() == *r1.read().unwrap()
            }
            (Neg(l0), Neg(r0)) => *l0.read().unwrap() == *r0.read().unwrap(),
            (Pow(l0, l1), Pow(r0, r1)) => *l0.read().unwrap() == *r0.read().unwrap() && l1 == r1,
            (Leaf(l0), Leaf(r0)) => l0 == r0,
            _ => false,
        }
    }
}

impl<T: Ord> PartialOrd for Expr<T> {
    fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<T: Ord> Ord for Expr<T> {
    fn cmp(&self, other: &Self) -> core::cmp::Ordering {
        use Expr::*;
        use Ordering::*;
        match (self, other) {
            (Pow(a, b), Pow(c, d)) => a
                .read()
                .unwrap()
                .cmp(&c.read().unwrap())
                .then_with(|| b.cmp(d)),
            (Leaf(a), Leaf(b)) => a.cmp(b),
            (Neg(a), Neg(b)) => a.read().unwrap().cmp(&b.read().unwrap()),
            (Add(a, b), Add(c, d)) | (Mul(a, b), Mul(c, d)) | (Div(a, b), Div(c, d)) => a
                .read()
                .unwrap()
                .cmp(&c.read().unwrap())
                .then_with(|| b.read().unwrap().cmp(&d.read().unwrap())),
            (_, Leaf(_)) => Greater,
            (Leaf(_), _) => Less,
            (_, Add(_, _)) => Greater,
            (Add(_, _), _) => Less,
            (_, Neg(_)) => Greater,
            (Neg(_), _) => Less,
            (_, Mul(_, _)) => Greater,
            (Mul(_, _), _) => Less,
            (_, Div(_, _)) => Greater,
            (Div(_, _), _) => Less,
        }
    }
}

impl<T> From<T> for Expr<T> {
    fn from(value: T) -> Self {
        Self::Leaf(value)
    }
}

impl<T: Zero> Sum<Self> for Expr<T> {
    fn sum<I: Iterator<Item = Self>>(mut iter: I) -> Self {
        iter.next().map_or(Self::Leaf(Zero::zero()), |expr| {
            iter.fold(expr, |a, b| a + b)
        })
    }
}

impl<T: One> Product<Self> for Expr<T> {
    fn product<I: Iterator<Item = Self>>(mut iter: I) -> Self {
        // TODO: zero or one?
        iter.next()
            .map_or(Self::Leaf(One::one()), |expr| iter.fold(expr, |a, b| a * b))
    }
}

impl<T: Clone> Mul<&Expr<T>> for &Expr<T> {
    type Output = Expr<T>;

    fn mul(self, rhs: &Expr<T>) -> Self::Output {
        Mul::mul(&(*self).clone(), rhs.clone())
    }
}

impl<T> Mul<Self> for Expr<T> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        Self::Mul(Arc::new(RwLock::new(self)), Arc::new(RwLock::new(rhs)))
    }
}

impl<T: Clone> Div<&Expr<T>> for &Expr<T> {
    type Output = Expr<T>;

    fn div(self, rhs: &Expr<T>) -> Self::Output {
        Div::div(self.clone(), rhs.clone())
    }
}

impl<T> Div<Self> for Expr<T> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self {
        Self::Div(Arc::new(RwLock::new(self)), Arc::new(RwLock::new(rhs)))
    }
}

impl<T: Clone> Add<&Expr<T>> for &Expr<T> {
    type Output = Expr<T>;

    fn add(self, rhs: &Expr<T>) -> Self::Output {
        Add::add(self.clone(), rhs.clone())
    }
}

impl<T> Add<Self> for Expr<T> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self::Add(Arc::new(RwLock::new(self)), Arc::new(RwLock::new(rhs)))
    }
}

impl<T: Clone> Sub<&Expr<T>> for &Expr<T> {
    type Output = Expr<T>;

    fn sub(self, rhs: &Expr<T>) -> Self::Output {
        Sub::sub(self.clone(), rhs.clone())
    }
}

impl<T> Sub<Self> for Expr<T> {
    type Output = Self;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn sub(self, rhs: Self) -> Self {
        self + rhs.neg()
    }
}

impl<T> Neg for Expr<T> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::Neg(Arc::new(RwLock::new(self)))
    }
}

impl<T: Clone> Neg for &Expr<T> {
    type Output = Expr<T>;

    #[inline]
    fn neg(self) -> Self::Output {
        self.clone().neg()
    }
}

impl<T> Mul<T> for Expr<T> {
    type Output = Self;

    fn mul(self, rhs: T) -> Self::Output {
        self * Self::Leaf(rhs)
    }
}

impl<T: Clone> Mul<&T> for &Expr<T> {
    type Output = Expr<T>;

    fn mul(self, rhs: &T) -> Self::Output {
        self * Expr::Leaf(rhs.clone())
    }
}

impl<T> Div<T> for Expr<T> {
    type Output = Self;

    fn div(self, rhs: T) -> Self::Output {
        self / Self::Leaf(rhs)
    }
}

impl<T: Clone> Div<&T> for &Expr<T> {
    type Output = Expr<T>;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn div(self, rhs: &T) -> Self::Output {
        self / Expr::Leaf(rhs.clone())
    }
}

impl<T> Add<T> for Expr<T> {
    type Output = Self;

    fn add(self, rhs: T) -> Self::Output {
        self + Self::Leaf(rhs)
    }
}

impl<T: Clone> Add<&T> for &Expr<T> {
    type Output = Expr<T>;

    fn add(self, rhs: &T) -> Self::Output {
        self + Expr::Leaf(rhs.clone())
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl<T> Sub<T> for Expr<T> {
    type Output = Self;

    fn sub(self, rhs: T) -> Self::Output {
        self + Self::Neg(Arc::new(RwLock::new(Self::Leaf(rhs))))
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl<T: Clone> Sub<&T> for &Expr<T> {
    type Output = Expr<T>;

    fn sub(self, rhs: &T) -> Self::Output {
        self + Expr::Neg(Arc::new(RwLock::new(Expr::Leaf(rhs.clone()))))
    }
}

forward_ref_binop!(impl< T: Clone > Mul, mul for Expr<T>, Expr<T>);
forward_ref_binop!(impl< T: Clone > Div, div for Expr<T>, Expr<T>);
forward_ref_binop!(impl< T: Clone > Add, add for Expr<T>, Expr<T>);
forward_ref_binop!(impl< T: Clone > Sub, sub for Expr<T>, Expr<T>);
forward_ref_binop!(impl< T: Clone > Mul, mul for Expr<T>, T);
forward_ref_binop!(impl< T: Clone > Div, div for Expr<T>, T);
forward_ref_binop!(impl< T: Clone > Add, add for Expr<T>, T);
forward_ref_binop!(impl< T: Clone > Sub, sub for Expr<T>, T);

impl<T> Pow<usize> for Expr<T> {
    type Output = Self;

    fn pow(self, rhs: usize) -> Self::Output {
        Self::Pow(Arc::new(RwLock::new(self)), rhs)
    }
}

impl<T: Clone> Pow<usize> for &Expr<T> {
    type Output = Expr<T>;

    fn pow(self, rhs: usize) -> Self::Output {
        Expr::Pow(Arc::new(RwLock::new(self.clone())), rhs)
    }
}

// TODO: remove reliance on clone
impl<T: Clone> MulAssign<Self> for Expr<T> {
    fn mul_assign(&mut self, other: Self) {
        *self = &*self * other;
    }
}

impl<T: Clone> MulAssign<T> for Expr<T> {
    fn mul_assign(&mut self, rhs: T) {
        *self = &*self * rhs;
    }
}

impl<T: Clone> DivAssign<Self> for Expr<T> {
    fn div_assign(&mut self, other: Self) {
        *self = &*self / other;
    }
}

impl<T: Clone> DivAssign<T> for Expr<T> {
    fn div_assign(&mut self, rhs: T) {
        *self = &*self / rhs;
    }
}

impl<T: Clone> AddAssign<Self> for Expr<T> {
    fn add_assign(&mut self, other: Self) {
        *self = &*self + other;
    }
}

impl<T: Clone> AddAssign<T> for Expr<T> {
    fn add_assign(&mut self, rhs: T) {
        *self = &*self + rhs;
    }
}

impl<T: Clone> SubAssign<Self> for Expr<T> {
    fn sub_assign(&mut self, other: Self) {
        *self = &*self - other;
    }
}

impl<T: Clone> SubAssign<T> for Expr<T> {
    fn sub_assign(&mut self, rhs: T) {
        *self = &*self - rhs;
    }
}

// TODO: why does removing clone here work? (has warnings)
forward_ref_op_assign!(impl< T: Clone > MulAssign, mul_assign for Expr<T>, Expr<T>);
forward_ref_op_assign!(impl< T: Clone > DivAssign, div_assign for Expr<T>, Expr<T>);
forward_ref_op_assign!(impl< T: Clone > AddAssign, add_assign for Expr<T>, Expr<T>);
forward_ref_op_assign!(impl< T: Clone > SubAssign, sub_assign for Expr<T>, Expr<T>);
forward_ref_op_assign!(impl< T: Clone > MulAssign, mul_assign for Expr<T>, T);
forward_ref_op_assign!(impl< T: Clone > DivAssign, div_assign for Expr<T>, T);
forward_ref_op_assign!(impl< T: Clone > AddAssign, add_assign for Expr<T>, T);
forward_ref_op_assign!(impl< T: Clone > SubAssign, sub_assign for Expr<T>, T);

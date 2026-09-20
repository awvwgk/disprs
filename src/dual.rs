use std::ops::{Add, AddAssign, Div, Mul, Neg, Sub};

#[derive(Clone)]
pub struct Dual2 {
    pub value: f64,
    pub gradient: Vec<f64>,
    pub hessian: Vec<f64>,
}

impl Dual2 {
    pub fn constant(value: f64, size: usize) -> Self {
        Self {
            value,
            gradient: vec![0.0; size],
            hessian: vec![0.0; size * size],
        }
    }

    pub fn variable(value: f64, size: usize, index: usize) -> Self {
        let mut result = Self::constant(value, size);
        result.gradient[index] = 1.0;
        result
    }

    fn unary(mut self, value: f64, first: f64, second: f64) -> Self {
        let size = self.gradient.len();
        for (index, item) in self.hessian.iter_mut().enumerate() {
            let row = index / size;
            let column = index % size;
            *item = first * *item + second * self.gradient[row] * self.gradient[column];
        }
        self.gradient.iter_mut().for_each(|item| *item *= first);
        self.value = value;
        self
    }

    pub fn exp(self) -> Self {
        let value = self.value.exp();
        self.unary(value, value, value)
    }

    pub fn powf(self, exponent: f64) -> Self {
        let input = self.value;
        let value = input.powf(exponent);
        self.unary(
            value,
            exponent * input.powf(exponent - 1.0),
            exponent * (exponent - 1.0) * input.powf(exponent - 2.0),
        )
    }

    pub fn sqrt(self) -> Self {
        self.powf(0.5)
    }

    fn scale(mut self, factor: f64) -> Self {
        self.value *= factor;
        self.gradient.iter_mut().for_each(|value| *value *= factor);
        self.hessian.iter_mut().for_each(|value| *value *= factor);
        self
    }

    fn product(mut self, rhs: Self) -> Self {
        let size = self.gradient.len();
        let left_value = self.value;
        for (index, item) in self.hessian.iter_mut().enumerate() {
            let row = index / size;
            let column = index % size;
            *item = *item * rhs.value
                + self.gradient[row] * rhs.gradient[column]
                + self.gradient[column] * rhs.gradient[row]
                + left_value * rhs.hessian[index];
        }
        for (left, right) in self.gradient.iter_mut().zip(rhs.gradient) {
            *left = *left * rhs.value + left_value * right;
        }
        self.value *= rhs.value;
        self
    }
}

impl Add for Dual2 {
    type Output = Self;

    fn add(mut self, rhs: Self) -> Self {
        self += rhs;
        self
    }
}

impl AddAssign for Dual2 {
    fn add_assign(&mut self, rhs: Self) {
        self.value += rhs.value;
        for (left, right) in self.gradient.iter_mut().zip(rhs.gradient) {
            *left += right;
        }
        for (left, right) in self.hessian.iter_mut().zip(rhs.hessian) {
            *left += right;
        }
    }
}

impl Add<f64> for Dual2 {
    type Output = Self;

    fn add(mut self, rhs: f64) -> Self {
        self.value += rhs;
        self
    }
}

impl Sub for Dual2 {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        self + -rhs
    }
}

impl Sub<f64> for Dual2 {
    type Output = Self;

    fn sub(mut self, rhs: f64) -> Self {
        self.value -= rhs;
        self
    }
}

impl Neg for Dual2 {
    type Output = Self;

    fn neg(self) -> Self {
        self.scale(-1.0)
    }
}

impl Mul<f64> for Dual2 {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self {
        self.scale(rhs)
    }
}

impl Mul for Dual2 {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        self.product(rhs)
    }
}

impl Div for Dual2 {
    type Output = Self;

    fn div(self, rhs: Self) -> Self {
        self.product(rhs.powf(-1.0))
    }
}

impl Div<f64> for Dual2 {
    type Output = Self;

    fn div(self, rhs: f64) -> Self {
        self.scale(rhs.recip())
    }
}

impl Div<Dual2> for f64 {
    type Output = Dual2;

    fn div(self, rhs: Dual2) -> Dual2 {
        rhs.powf(-1.0) * self
    }
}

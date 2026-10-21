//! Minimal B-spline bases and NURBS evaluation in arbitrary dimension.
//!
//! The free functions assume valid input shapes; [`nurbs::NURBSCurve`] and
//! [`nurbs::NURBSSurface`] validate their data at construction and their parameters
//! at evaluation. The library has no dependencies; Kuva is used only by an example.
//! All free evaluation functions build on recursive [`bspline_basis`]. The
//! structs' `eval` methods use a private iterative Cox–de Boor implementation.
//!
//! With `n` control points and degree `p`, a knot vector has `n + p + 1`
//! entries. The structs evaluate on `[knots[p], knots[n]]`; raw bases retain
//! their full knot support outside that domain. Interior
//! knots are evaluated from the right, and the domain's right endpoint from the
//! left. In particular, clamped splines interpolate their endpoint control points
//! when the corresponding weights are positive.
//!
//! ```
//! use bsfun::nurbs::NURBSCurve;
//!
//! let curve = NURBSCurve::new(
//!     1,
//!     vec![0.0, 0.0, 1.0, 1.0],
//!     vec![1.0, 1.0],
//!     vec![vec![0.0, 2.0], vec![4.0, 6.0]],
//! )?;
//! assert_eq!(curve.eval(0.5)?, vec![2.0, 4.0]);
//! assert_eq!(curve.eval(1.0)?, vec![4.0, 6.0]);
//! # Ok::<(), &'static str>(())
//! ```

pub mod nurbs;

/// Evaluate B-spline basis `i` by direct recursive Cox–de Boor recurrence.
///
/// Returns zero outside its knot support, including for non-finite `t`.
/// Repeated knots contribute zero when their recurrence denominator is zero.
/// The standard domain's right endpoint uses its left-hand limit.
/// All free evaluation functions build on this recursive definition. The structs'
/// `eval` methods use a private iterative implementation to avoid repeated work.
///
/// # Panics
///
/// Panics if the degree or index cannot fit in the knot vector. Knots must be
/// finite, non-decreasing, and describe a nonempty domain, as required by
/// [`nurbs::NURBSCurve::new`]. The free functions do not validate these conditions.
pub fn bspline_basis(i: usize, degree: usize, knots: &[f64], t: f64) -> f64 {
    // Preserve the original endpoint convention throughout the recursion: the
    // standard domain depends on the requested degree, not each lower degree.
    fn recurse(i: usize, degree: usize, knots: &[f64], t: f64, from_left: bool) -> f64 {
        let in_support = if from_left {
            knots[i] < t && t <= knots[i + degree + 1]
        } else {
            knots[i] <= t && t < knots[i + degree + 1]
        };
        if !in_support {
            return 0.0;
        }
        if degree == 0 {
            return 1.0;
        }

        let left_basis = recurse(i, degree - 1, knots, t, from_left);
        let right_basis = recurse(i + 1, degree - 1, knots, t, from_left);
        let left = if left_basis == 0.0 || knots[i] == knots[i + degree] {
            0.0
        } else {
            knot_fraction(t, knots[i], knots[i + degree]) * left_basis
        };
        let right = if right_basis == 0.0 || knots[i + 1] == knots[i + degree + 1] {
            0.0
        } else {
            knot_fraction(-t, -knots[i + degree + 1], -knots[i + 1]) * right_basis
        };
        left + right
    }

    assert!(
        degree < knots.len().saturating_sub(1),
        "Invalid degree or knot count"
    );
    let count = knots.len() - degree - 1;
    assert!(i < count, "Invalid basis index");
    if !t.is_finite() {
        return 0.0;
    }
    let from_left = t == knots[count] || t == knots[knots.len() - 1];
    recurse(i, degree, knots, t, from_left)
}

/// Evaluate rational curve basis `i`, with each weight applied once.
///
/// Weights must be finite and non-negative. Returns zero outside knot support or
/// where the rational denominator is zero. Weight scaling avoids dependence on
/// a common weight magnitude (there is no absolute epsilon cutoff).
/// Polynomial bases are evaluated with recursive [`bspline_basis`].
///
/// # Panics
///
/// Input shapes and indices must match [`nurbs::NURBSCurve::new`]. Invalid shapes
/// may panic. Use the struct to validate data and detect undefined evaluations.
pub fn rational_bspline_basis_curve(
    i: usize,
    degree: usize,
    knots: &[f64],
    weights: &[f64],
    t: f64,
) -> f64 {
    recursive_curve_basis(degree, knots, weights, t)[i]
}

/// Evaluate rational tensor-product surface basis `(i, j)`.
///
/// Each weight multiplies the product of the two B-spline bases once. Returns
/// zero outside knot support or where the rational denominator is zero.
/// Both polynomial bases are evaluated with recursive [`bspline_basis`].
///
/// # Panics
///
/// Input shapes, values, and indices must match [`nurbs::NURBSSurface::new`].
/// Invalid shapes may panic. Use the struct for validated evaluation.
#[allow(clippy::too_many_arguments)] // Preserve the small, explicit public API.
pub fn rational_bspline_basis_surface(
    i: usize,
    j: usize,
    degree_u: usize,
    degree_v: usize,
    knots_u: &[f64],
    knots_v: &[f64],
    weights: &[Vec<f64>],
    u: f64,
    v: f64,
) -> f64 {
    let basis = recursive_surface_basis(degree_u, degree_v, knots_u, knots_v, weights, u, v);
    assert!(
        i < weights.len() && j < weights[0].len(),
        "Invalid basis index"
    );
    basis[i * weights[0].len() + j]
}

/// Evaluate a NURBS curve point in any positive dimension.
///
/// Returns a zero vector outside knot support or where all supporting weights
/// are zero. Unlike [`nurbs::NURBSCurve::eval`], this function does not restrict
/// parameters to the standard evaluation domain.
/// Uses recursive [`bspline_basis`] and computes the rational denominator once.
///
/// # Panics
///
/// Input shapes and values must match [`nurbs::NURBSCurve::new`]. Invalid shapes
/// may panic. Data validation is deliberately left to the struct constructor.
pub fn nurbs_curve_point(
    t: f64,
    control_points: &[Vec<f64>],
    weights: &[f64],
    knots: &[f64],
    degree: usize,
) -> Vec<f64> {
    let basis = recursive_curve_basis(degree, knots, weights, t);
    weighted_point(&basis, control_points.iter().map(Vec::as_slice))
}

/// Evaluate a tensor-product NURBS surface point in any positive dimension.
///
/// Returns a zero vector outside knot support or where all supporting weights
/// are zero. Unlike [`nurbs::NURBSSurface::eval`], this function does not restrict
/// parameters to the standard evaluation domains.
/// Uses recursive [`bspline_basis`] and computes the rational denominator once.
///
/// # Panics
///
/// Input shapes and values must match [`nurbs::NURBSSurface::new`]. Invalid shapes
/// may panic. Data validation is deliberately left to the struct constructor.
#[allow(clippy::too_many_arguments)]
pub fn nurbs_surface_point(
    control_points: &[Vec<Vec<f64>>],
    weights: &[Vec<f64>],
    knots_u: &[f64],
    knots_v: &[f64],
    degree_u: usize,
    degree_v: usize,
    u: f64,
    v: f64,
) -> Vec<f64> {
    let basis = recursive_surface_basis(degree_u, degree_v, knots_u, knots_v, weights, u, v);
    weighted_point(&basis, control_points.iter().flatten().map(Vec::as_slice))
}

// The free functions follow the recursive basis definition. Keep normalization
// and point summation shared with the iterative path for the same numerics.
fn recursive_curve_basis(degree: usize, knots: &[f64], weights: &[f64], t: f64) -> Vec<f64> {
    let mut basis = weights
        .iter()
        .enumerate()
        .map(|(i, _)| bspline_basis(i, degree, knots, t))
        .collect::<Vec<_>>();
    rationalize(&mut basis, weights.iter().copied());
    basis
}

#[allow(clippy::too_many_arguments)]
fn recursive_surface_basis(
    degree_u: usize,
    degree_v: usize,
    knots_u: &[f64],
    knots_v: &[f64],
    weights: &[Vec<f64>],
    u: f64,
    v: f64,
) -> Vec<f64> {
    let columns = weights.first().map_or(0, Vec::len);
    let basis_v: Vec<_> = (0..columns)
        .map(|j| bspline_basis(j, degree_v, knots_v, v))
        .collect();
    let mut basis: Vec<_> = weights
        .iter()
        .enumerate()
        .flat_map(|(i, row)| {
            assert_eq!(row.len(), columns, "Invalid weights columns");
            let basis_u = bspline_basis(i, degree_u, knots_u, u);
            basis_v.iter().map(move |basis_v| basis_u * basis_v)
        })
        .collect();
    rationalize(&mut basis, weights.iter().flatten().copied());
    basis
}

// The structs evaluate all bases together, bottom up: O(p * (n + p)), without recursion or
// repeatedly evaluating the same lower-degree bases. Updating left to right is
// safe because each entry only reads itself and its still-unmodified neighbour.
fn basis_values(degree: usize, knots: &[f64], t: f64) -> Vec<f64> {
    assert!(
        degree < knots.len().saturating_sub(1),
        "Invalid degree or knot count"
    );
    let count = knots.len() - degree - 1;
    let mut basis = vec![0.0; knots.len() - 1];
    if t.is_finite() && t >= knots[0] && t <= knots[knots.len() - 1] {
        // At the right endpoint choose the last nonempty span to its left,
        // including when that endpoint has repeated knots or is not clamped.
        let span = if t == knots[count] || t == knots[knots.len() - 1] {
            knots.partition_point(|&knot| knot < t) - 1
        } else {
            knots.partition_point(|&knot| knot <= t) - 1
        };
        basis[span] = 1.0;
        for order in 1..=degree {
            for i in 0..knots.len() - order - 1 {
                let left = if basis[i] == 0.0 || knots[i] == knots[i + order] {
                    0.0
                } else {
                    knot_fraction(t, knots[i], knots[i + order]) * basis[i]
                };
                let right = if basis[i + 1] == 0.0 || knots[i + 1] == knots[i + order + 1] {
                    0.0
                } else {
                    // Reflect the interval to compute (b - t) / (b - a)
                    // directly, avoiding cancellation from 1 - fraction.
                    knot_fraction(-t, -knots[i + order + 1], -knots[i + 1]) * basis[i + 1]
                };
                basis[i] = left + right;
            }
        }
    }
    basis.truncate(count);
    basis
}

// (x - a) / (b - a), including finite endpoints whose difference overflows.
fn knot_fraction(x: f64, a: f64, b: f64) -> f64 {
    let width = b - a;
    if width.is_infinite() {
        (x * 0.5 - a * 0.5) / (b * 0.5 - a * 0.5)
    } else {
        (x - a) / width
    }
}

fn curve_basis(degree: usize, knots: &[f64], weights: &[f64], t: f64) -> Vec<f64> {
    let mut basis = basis_values(degree, knots, t);
    assert_eq!(basis.len(), weights.len(), "Invalid weights length");
    rationalize(&mut basis, weights.iter().copied());
    basis
}

#[allow(clippy::too_many_arguments)]
fn surface_basis(
    degree_u: usize,
    degree_v: usize,
    knots_u: &[f64],
    knots_v: &[f64],
    weights: &[Vec<f64>],
    u: f64,
    v: f64,
) -> Vec<f64> {
    let basis_u = basis_values(degree_u, knots_u, u);
    let basis_v = basis_values(degree_v, knots_v, v);
    assert_eq!(weights.len(), basis_u.len(), "Invalid weights rows");
    assert!(
        weights.iter().all(|row| row.len() == basis_v.len()),
        "Invalid weights columns"
    );
    let mut basis: Vec<_> = basis_u
        .iter()
        .flat_map(|bu| basis_v.iter().map(move |bv| bu * bv))
        .collect();
    rationalize(&mut basis, weights.iter().flatten().copied());
    basis
}

// Scale the weighted contributions, not just the weights: a large weight can
// have a tiny basis, and vice versa. Their products may still be comparable.
fn rationalize(basis: &mut [f64], weights: impl Iterator<Item = f64> + Clone) {
    let scale = basis
        .iter()
        .zip(weights.clone())
        .map(|(value, weight)| value * weight)
        .fold(0.0_f64, f64::max);
    if scale == 0.0 {
        // All products may have underflowed. Rescale supporting weights first
        // so at least one positive contribution survives, if there is one.
        let weight_scale = basis
            .iter()
            .zip(weights.clone())
            .filter(|(value, _)| **value > 0.0)
            .map(|(_, weight)| weight)
            .fold(0.0_f64, f64::max);
        if weight_scale == 0.0 {
            basis.fill(0.0);
            return;
        }
        for (value, weight) in basis.iter_mut().zip(weights) {
            if *value != 0.0 {
                *value *= weight / weight_scale;
            }
        }
    } else {
        for (value, weight) in basis.iter_mut().zip(weights) {
            *value = if scale >= f64::MIN_POSITIVE {
                (*value * weight) / scale
            } else if *value < weight {
                // For subnormal products, divide the smaller factor first to
                // retain precision without overflowing the intermediate ratio.
                (*value / scale) * weight
            } else {
                (weight / scale) * *value
            };
        }
    }
    let denominator = compensated_sum(basis.iter().copied());
    if denominator != 0.0 {
        for value in basis {
            *value /= denominator;
        }
    }
}

fn weighted_point<'a>(basis: &[f64], points: impl Iterator<Item = &'a [f64]> + Clone) -> Vec<f64> {
    let dimensions = points
        .clone()
        .next()
        .expect("Control points must not be empty")
        .len();
    assert_eq!(
        basis.len(),
        points.clone().count(),
        "Invalid control point count"
    );
    (0..dimensions)
        .map(|axis| {
            let scale = points
                .clone()
                .zip(basis)
                .filter(|(_, value)| **value > 0.0)
                .map(|(point, _)| point[axis].abs())
                .fold(0.0_f64, f64::max);
            if scale == 0.0 {
                return 0.0;
            }
            let coordinate = compensated_sum(
                points
                    .clone()
                    .zip(basis)
                    .filter(|(_, value)| **value > 0.0)
                    .map(|(point, value)| (point[axis] / scale) * value),
            );
            // A convex combination stays within [-scale, scale]. Clamp tiny
            // roundoff excursions before rescaling, which could overflow at MAX.
            coordinate.clamp(-1.0, 1.0) * scale
        })
        .collect()
}

// Neumaier summation retains small terms that would otherwise be lost between
// large opposite contributions. Inputs here are scaled to avoid overflow.
fn compensated_sum(values: impl Iterator<Item = f64>) -> f64 {
    let (sum, correction) = values.fold((0.0_f64, 0.0), |(sum, correction), value| {
        let next = sum + value;
        let error = if sum.abs() >= value.abs() {
            (sum - next) + value
        } else {
            (value - next) + sum
        };
        (next, correction + error)
    });
    sum + correction
}

fn parameter_in_domain(t: f64, degree: usize, knots: &[f64]) -> bool {
    t.is_finite() && t >= knots[degree] && t <= knots[knots.len() - degree - 1]
}

fn validate_knots(degree: usize, knots: &[f64], count: usize) -> Result<(), &'static str> {
    if degree >= count {
        return Err("Degree must be less than the number of control points");
    }
    if count.checked_add(degree).and_then(|n| n.checked_add(1)) != Some(knots.len()) {
        return Err("Number of knots must be number of control points + degree + 1");
    }
    if !knots.iter().all(|knot| knot.is_finite()) {
        return Err("Knots must be finite");
    }
    if knots.windows(2).any(|pair| pair[0] > pair[1]) {
        return Err("Knots must be non-decreasing");
    }
    if knots[degree] >= knots[count] {
        return Err("The evaluation domain must have positive length");
    }
    let mut multiplicity = 1;
    for pair in knots.windows(2) {
        multiplicity = if pair[0] == pair[1] {
            multiplicity + 1
        } else {
            1
        };
        if multiplicity > degree + 1 {
            return Err("Knot multiplicity must not exceed degree + 1");
        }
    }
    Ok(())
}

fn validate_weights(weights: impl Iterator<Item = f64>) -> Result<(), &'static str> {
    let mut any_positive = false;
    for weight in weights {
        if !weight.is_finite() || weight < 0.0 {
            return Err("Weights must be finite and non-negative");
        }
        any_positive |= weight > 0.0;
    }
    if !any_positive {
        return Err("At least one weight must be positive");
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{basis_values, bspline_basis};

    #[test]
    fn recursive_and_iterative_bases_agree_on_support_and_endpoint_conventions() {
        let mut cases = vec![(2, vec![0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 5.0])];
        for degree in 0..=5 {
            let clamped = vec![-2.0; degree + 1]
                .into_iter()
                .chain([-1.0, 0.25, 0.5])
                .chain(vec![2.0; degree + 1])
                .collect();
            let unclamped = (0..2 * degree + 5).map(|i| i as f64 * 0.5).collect();
            cases.extend([(degree, clamped), (degree, unclamped)]);
        }
        for (degree, knots) in cases {
            let parameters = knots
                .iter()
                .copied()
                .chain(knots.windows(2).map(|pair| (pair[0] + pair[1]) / 2.0))
                .chain([-3.0, 20.0, f64::NAN, f64::NEG_INFINITY, f64::INFINITY]);
            for t in parameters {
                for (i, expected) in basis_values(degree, &knots, t).into_iter().enumerate() {
                    let actual = bspline_basis(i, degree, &knots, t);
                    assert!(
                        (actual - expected).abs() < 1e-14,
                        "basis {i}, degree {degree}, t={t}: {actual} != {expected}"
                    );
                }
            }
        }
    }
}

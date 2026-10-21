use crate::{curve_basis, parameter_in_domain, validate_knots, validate_weights, weighted_point};

/// A validated NURBS curve in arbitrary positive dimension.
///
/// Stores its input data unchanged. Evaluation uses a private iterative
/// Cox–de Boor implementation; the free functions use recursive bases. Both
/// paths share weight normalization and coordinate summation.
pub struct NURBSCurve {
    degree: usize,
    knots: Vec<f64>,
    weights: Vec<f64>,
    control_points: Vec<Vec<f64>>,
}

impl NURBSCurve {
    /// Construct a curve, validating knot counts/order/multiplicity, a nonempty
    /// domain, finite coordinates of equal positive dimension, and finite
    /// non-negative weights with at least one positive weight.
    ///
    /// Degree zero is supported; the degree must be less than the point count.
    /// Zero weights are allowed, but can leave some parameters undefined.
    pub fn new(
        degree: usize,
        knots: Vec<f64>,
        weights: Vec<f64>,
        control_points: Vec<Vec<f64>>,
    ) -> Result<Self, &'static str> {
        let dimensions = control_points.first().map_or(0, Vec::len);
        if dimensions == 0 {
            return Err("Control points must be nonempty and have positive dimension");
        }
        if control_points.iter().any(|point| {
            point.len() != dimensions || point.iter().any(|coordinate| !coordinate.is_finite())
        }) {
            return Err("Control points must have equal dimension and finite coordinates");
        }
        if weights.len() != control_points.len() {
            return Err("Number of weights and control points must be the same");
        }
        validate_knots(degree, &knots, control_points.len())?;
        validate_weights(weights.iter().copied())?;
        Ok(Self {
            degree,
            knots,
            weights,
            control_points,
        })
    }

    /// Evaluate on `[knots[degree], knots[control_point_count]]`.
    /// Uses the private iterative Cox–de Boor implementation.
    ///
    /// Returns an error for non-finite or out-of-domain `t`, or a zero rational
    /// denominator. The right endpoint is evaluated by its left-hand limit.
    pub fn eval(&self, t: f64) -> Result<Vec<f64>, &'static str> {
        if !parameter_in_domain(t, self.degree, &self.knots) {
            return Err("Parameter t must be finite and within the evaluation domain");
        }
        let basis = curve_basis(self.degree, &self.knots, &self.weights, t);
        if !basis.iter().any(|&value| value > 0.0) {
            return Err("The rational denominator is zero at this parameter");
        }
        Ok(weighted_point(
            &basis,
            self.control_points.iter().map(Vec::as_slice),
        ))
    }
}

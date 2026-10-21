use crate::{parameter_in_domain, surface_basis, validate_knots, validate_weights, weighted_point};

/// A validated tensor-product NURBS surface in arbitrary positive dimension.
///
/// Rows follow the `u` direction and columns follow `v`, for both weights and
/// control points. All rows must have the same length.
/// Evaluation uses a private iterative Cox–de Boor implementation; the free
/// functions use recursive bases. Both paths share weight normalization and
/// coordinate summation.
pub struct NURBSSurface {
    degree_u: usize,
    degree_v: usize,
    knots_u: Vec<f64>,
    knots_v: Vec<f64>,
    weights: Vec<Vec<f64>>,
    control_points: Vec<Vec<Vec<f64>>>,
}

impl NURBSSurface {
    /// Construct a surface, validating a rectangular control net and weight
    /// grid, equal positive point dimensions, finite coordinates and knots,
    /// non-negative finite weights (at least one positive), and valid knot
    /// counts, order, multiplicities, and nonempty domains in both directions.
    ///
    /// Degree zero is supported in either direction. Each degree must be less
    /// than the corresponding control point count. Zero weights can leave some
    /// parameter pairs undefined.
    pub fn new(
        degree_u: usize,
        degree_v: usize,
        knots_u: Vec<f64>,
        knots_v: Vec<f64>,
        weights: Vec<Vec<f64>>,
        control_points: Vec<Vec<Vec<f64>>>,
    ) -> Result<Self, &'static str> {
        let rows = control_points.len();
        let columns = control_points.first().map_or(0, Vec::len);
        let dimensions = control_points
            .first()
            .and_then(|row| row.first())
            .map_or(0, Vec::len);
        if dimensions == 0 {
            return Err("The control net must be nonempty with positive point dimension");
        }
        if control_points.iter().any(|row| row.len() != columns) {
            return Err("The control net must be rectangular");
        }
        if control_points.iter().flatten().any(|point| {
            point.len() != dimensions || point.iter().any(|coordinate| !coordinate.is_finite())
        }) {
            return Err("Control points must have equal dimension and finite coordinates");
        }
        if weights.len() != rows || weights.iter().any(|row| row.len() != columns) {
            return Err("Weights must match the rectangular control net");
        }
        validate_knots(degree_u, &knots_u, rows)?;
        validate_knots(degree_v, &knots_v, columns)?;
        validate_weights(weights.iter().flatten().copied())?;
        Ok(Self {
            degree_u,
            degree_v,
            knots_u,
            knots_v,
            weights,
            control_points,
        })
    }

    /// Evaluate on the closed knot domains in `u` and `v`.
    /// Uses the private iterative Cox–de Boor implementation.
    ///
    /// Each domain runs from `knots[degree]` to `knots[control_point_count]`.
    /// Returns an error for non-finite or out-of-domain parameters, or a zero
    /// rational denominator. Right endpoints use their left-hand limits.
    pub fn eval(&self, u: f64, v: f64) -> Result<Vec<f64>, &'static str> {
        if !parameter_in_domain(u, self.degree_u, &self.knots_u) {
            return Err("Parameter u must be finite and within the evaluation domain");
        }
        if !parameter_in_domain(v, self.degree_v, &self.knots_v) {
            return Err("Parameter v must be finite and within the evaluation domain");
        }
        let basis = surface_basis(
            self.degree_u,
            self.degree_v,
            &self.knots_u,
            &self.knots_v,
            &self.weights,
            u,
            v,
        );
        if !basis.iter().any(|&value| value > 0.0) {
            return Err("The rational denominator is zero at these parameters");
        }
        Ok(weighted_point(
            &basis,
            self.control_points.iter().flatten().map(Vec::as_slice),
        ))
    }
}

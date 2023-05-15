use crate::rational_bspline_basis_curve;

pub struct NURBSCurve {
    degree: usize,
    knots: Vec<f64>,
    weights: Vec<f64>,
    control_points: Vec<Vec<f64>>,
}

impl NURBSCurve {
    pub fn new(
        degree: usize,
        knots: Vec<f64>,
        weights: Vec<f64>,
        control_points: Vec<Vec<f64>>,
    ) -> Result<Self, &'static str> {
        // Check the validity of the inputs
        if degree == 0 {
            return Err("Degree must be greater than 0");
        }
        if weights.len() != control_points.len() {
            return Err("Number of weights and control points must be the same");
        }
        if knots.len() != control_points.len() + degree + 1 {
            return Err("Number of knots must be number of control points + degree + 1");
        }
        // Check that the knots are non-decreasing
        for i in 0..knots.len() - 1 {
            if knots[i] > knots[i + 1] {
                return Err("Knots must be non-decreasing");
            }
        }
        // Check that the weights are non-negative
        for &w in &weights {
            if w < 0.0 {
                return Err("Weights must be non-negative");
            }
        }
        // the implementation uses Vec<Vec<>> so points could potentially be of different dimension
        if !control_points.iter().all(|ref v| v.len() == control_points[0].len()) {
            return Err("All control points must be of the same dimension");
        }

        // Everything checks out, so construct the NURBSCurve
        Ok(NURBSCurve {
            degree,
            knots,
            weights,
            control_points,
        })
    }

    pub fn eval(&self, t: f64) -> Result<Vec<f64>, &'static str> {
        if t < *self.knots.first().unwrap_or(&0.0) || t > *self.knots.last().unwrap_or(&0.0) {
            return Err("The parameter t is out of range of the knot values");
        }
        let n = self.control_points.len();
        let d = self.control_points[0].len();
        let mut result = vec![0.0; d];
        for i in 0..n {
            let basis = rational_bspline_basis_curve(i, self.degree, &self.knots, &self.weights, t);
            for j in 0..d {
                // fill every axis of the result vector
                result[j] += basis * self.control_points[i][j];
            }
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nurbs_curve() {
        // Define the degree, knots, weights, and control points
        let degree = 2;
        let knots = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        let weights = vec![1.0, 1.0, 1.0];
        let control_points = vec![vec![0.0, 0.0], vec![0.5, 0.5], vec![1.0, 0.0]];

        // Construct the NURBSCurve
        let curve = NURBSCurve::new(degree, knots, weights, control_points).unwrap();

        // Evaluate the curve at a few parameter values
        let t_values = vec![0.0, 0.25, 0.5, 0.75, 1.0];
        for &t in &t_values {
            let point = curve.eval(t).unwrap();
            println!("Curve at t = {}: {:?}", t, point);
        }
    }
}

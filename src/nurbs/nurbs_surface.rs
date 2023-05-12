use crate::rational_bspline_basis_surface;
pub struct NURBSSurface {
    degree_u: usize,
    degree_v: usize,
    knots_u: Vec<f64>,
    knots_v: Vec<f64>,
    weights: Vec<Vec<f64>>,
    control_points: Vec<Vec<Vec<f64>>>, // each control point is a 3D point
}

impl NURBSSurface {
    pub fn new(degree_u: usize, degree_v: usize, knots_u: Vec<f64>, knots_v: Vec<f64>, weights: Vec<Vec<f64>>, control_points: Vec<Vec<Vec<f64>>>) -> Result<Self, &'static str> {
        if knots_u.len() != control_points.len() + degree_u + 1 {
            return Err("Invalid knots_u length");
        }
        if knots_v.len() != control_points[0].len() + degree_v + 1 {
            return Err("Invalid knots_v length");
        }
        if weights.len() != control_points.len() || weights[0].len() != control_points[0].len() {
            return Err("Invalid weights dimensions");
        }
        Ok(Self { degree_u, degree_v, knots_u, knots_v, weights, control_points })
    }

    pub fn eval(&self, u: f64, v: f64) -> Result<Vec<f64>, &'static str> {
        if u < *self.knots_u.first().unwrap_or(&0.0) || u > *self.knots_u.last().unwrap_or(&0.0) {
            return Err("Parameter u is out of the knot spans bounds");
        }
        if v < *self.knots_v.first().unwrap_or(&0.0) || v > *self.knots_v.last().unwrap_or(&0.0) {
            return Err("Parameter v is out of the knot spans bounds");
        }

        let n = self.control_points.len();
        let m = self.control_points[0].len();
        let d = self.control_points[0][0].len();

        let mut result = vec![0.0; d];

        for i in 0..n {
            for j in 0..m {
                let basis = rational_bspline_basis_surface(i, j, self.degree_u, self.degree_v, &self.knots_u, &self.knots_v, &self.weights, u, v);
                for k in 0..d {
                    result[k] += basis * self.control_points[i][j][k];
                }
            }
        }

        Ok(result)
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nurbs_surface() {
        let degree_u = 1;
        let degree_v = 1;

        let knots_u = vec![0.0, 0.0, 1.0, 1.0];
        let knots_v = vec![0.0, 0.0, 1.0, 1.0];

        let weights = vec![
            vec![1.0, 1.0],
            vec![1.0, 1.0],
        ];

        let control_points = vec![
            vec![
                vec![0.0, 0.0, 0.0], // Point (0, 0, 0)
                vec![1.0, 0.0, 0.0], // Point (1, 0, 0)
            ],
            vec![
                vec![0.0, 1.0, 0.0], // Point (0, 1, 0)
                vec![1.0, 1.0, 0.0], // Point (1, 1, 0)
            ],
        ];

        let surface = NURBSSurface::new(degree_u, degree_v, knots_u, knots_v, weights, control_points).unwrap();

        let u = 0.5;
        let v = 0.75;

        let result = surface.eval(u, v).unwrap();
        assert_eq!(result, vec![0.75, 0.5, 0.0]);  // Since we're on the xy plane, z should be 0
    }
}

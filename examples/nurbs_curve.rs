use bsfun::nurbs_curve_point;

fn main() {
    // Define some 3D control points.
    let control_points = vec![
        vec![0.0, 0.0, 0.0],
        vec![0.5, 0.5, 0.5],
        vec![1.0, 0.0, 0.0],
    ];

    // Define some weights.
    let weights = vec![1.0, 0.5, 1.0];

    // Define some knots.
    let knots = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];

    // Define the degree.
    let degree = 2;

    // Calculate points on the NURBS curve.
    for i in 0..=10 {
        let t = i as f64 / 10.0;
        let point = nurbs_curve_point(t, &control_points, &weights, &knots, degree);
        println!("Point at t = {}: {:?}", t, point);
    }
}

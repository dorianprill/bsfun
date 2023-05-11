use bsfun::nurbs_surface_point;

fn main() {
    // Define some 3D control points.
    let control_points = vec![
        vec![
            vec![0.0, 0.0, 0.0],
            vec![0.0, 0.5, 0.5],
            vec![0.0, 1.0, 0.0],
        ],
        vec![
            vec![0.5, 0.0, 0.5],
            vec![0.5, 0.5, 1.0],
            vec![0.5, 1.0, 0.5],
        ],
        vec![
            vec![1.0, 0.0, 0.0],
            vec![1.0, 0.5, 0.5],
            vec![1.0, 1.0, 0.0],
        ],
    ];

    // Define some weights.
    let weights = vec![
        vec![1.0, 0.5, 1.0],
        vec![0.5, 1.0, 0.5],
        vec![1.0, 0.5, 1.0],
    ];

    // Define some knots.
    let knots_u = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
    let knots_v = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];

    // Define the degree.
    let degree_u = 2;
    let degree_v = 2;

    // Calculate points on the NURBS surface.
    for i in 0..=10 {
        for j in 0..=10 {
            let u = i as f64 / 10.0;
            let v = j as f64 / 10.0;
            let point = nurbs_surface_point(
                u,
                v,
                &control_points,
                &weights,
                &knots_u,
                &knots_v,
                degree_u,
                degree_v,
            );
            println!("Point at (u, v) = ({}, {}): {:?}", u, v, point);
        }
    }
}

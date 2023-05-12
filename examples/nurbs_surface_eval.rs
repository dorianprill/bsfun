use bsfun::nurbs::NURBSSurface;

fn main() {

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

    println!("Point on surface: {:?} (Should be [0.75, 0.5, 0.0])", result);

}
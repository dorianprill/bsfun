use bsfun::nurbs::{NURBSCurve, NURBSSurface};

type Point = Vec<f64>;
type CurveInput = (usize, Vec<f64>, Vec<f64>, Vec<Point>);

fn curve_input() -> CurveInput {
    (
        1,
        vec![0.0, 0.0, 1.0, 1.0],
        vec![1.0, 1.0],
        vec![vec![0.0], vec![1.0]],
    )
}

#[test]
fn curve_constructor_rejects_invalid_shapes_without_panicking() {
    let invalid: Vec<CurveInput> = vec![
        (1, vec![], vec![], vec![]),
        (1, vec![0.0, 0.0, 1.0, 1.0], vec![], vec![]),
        (0, vec![0.0, 1.0], vec![1.0], vec![vec![]]),
        (
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![1.0],
            vec![vec![0.0], vec![1.0]],
        ),
        (
            1,
            vec![0.0, 1.0],
            vec![1.0, 1.0],
            vec![vec![0.0], vec![1.0]],
        ),
        (
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![1.0, 1.0],
            vec![vec![0.0], vec![1.0, 2.0]],
        ),
    ];
    for (degree, knots, weights, points) in invalid {
        assert!(NURBSCurve::new(degree, knots, weights, points).is_err());
    }
    for degree in [2, usize::MAX - 1, usize::MAX] {
        let (_, knots, weights, points) = curve_input();
        assert!(NURBSCurve::new(degree, knots, weights, points).is_err());
    }
}

#[test]
fn curve_constructor_rejects_invalid_knots_weights_and_coordinates() {
    for knots in [
        vec![0.0, 0.5, 0.25, 1.0],
        vec![1.0; 4],
        vec![0.0, 0.0, 0.0, 1.0],
    ] {
        let (degree, _, weights, points) = curve_input();
        assert!(NURBSCurve::new(degree, knots, weights, points).is_err());
    }
    // Excess multiplicity even with a positive evaluation domain.
    assert!(NURBSCurve::new(
        1,
        vec![0.0, 0.0, 0.0, 1.0, 1.0],
        vec![1.0; 3],
        vec![vec![0.0]; 3]
    )
    .is_err());
    for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        for i in 0..4 {
            let (degree, mut knots, weights, points) = curve_input();
            knots[i] = invalid;
            assert!(NURBSCurve::new(degree, knots, weights, points).is_err());
        }
        for i in 0..2 {
            let (degree, knots, mut weights, points) = curve_input();
            weights[i] = invalid;
            assert!(NURBSCurve::new(degree, knots, weights, points).is_err());
            let (degree, knots, weights, mut points) = curve_input();
            points[i][0] = invalid;
            assert!(NURBSCurve::new(degree, knots, weights, points).is_err());
        }
    }
    for weights in [vec![-1.0, 1.0], vec![0.0, 0.0]] {
        let (degree, knots, _, points) = curve_input();
        assert!(NURBSCurve::new(degree, knots, weights, points).is_err());
    }
}

fn surface_with_data(
    weights: Vec<Vec<f64>>,
    points: Vec<Vec<Point>>,
) -> Result<NURBSSurface, &'static str> {
    NURBSSurface::new(
        1,
        1,
        vec![0.0, 0.0, 1.0, 1.0],
        vec![0.0, 0.0, 1.0, 1.0],
        weights,
        points,
    )
}

#[test]
fn surface_constructor_rejects_empty_and_ragged_nets_and_weight_grids() {
    for points in [
        vec![],
        vec![vec![]],
        vec![vec![vec![]]],
        vec![vec![vec![0.0]], vec![]],
        vec![vec![vec![0.0]; 2], vec![vec![0.0]; 1]],
        vec![vec![vec![0.0]; 2], vec![vec![0.0]; 3]],
        vec![vec![vec![0.0]; 2], vec![vec![0.0], vec![0.0, 1.0]]],
    ] {
        assert!(surface_with_data(vec![vec![1.0; 2]; 2], points).is_err());
    }
    for weights in [
        vec![],
        vec![vec![]],
        vec![vec![1.0; 2]],
        vec![vec![1.0; 2], vec![]],
        vec![vec![1.0; 2], vec![1.0]],
        vec![vec![1.0; 2], vec![1.0; 3]],
    ] {
        assert!(surface_with_data(weights, vec![vec![vec![0.0]; 2]; 2]).is_err());
    }
}

#[test]
fn surface_constructor_checks_every_coordinate_and_weight() {
    for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        for i in 0..2 {
            for j in 0..2 {
                for axis in 0..3 {
                    let mut points = vec![vec![vec![0.0; 3]; 2]; 2];
                    points[i][j][axis] = invalid;
                    assert!(surface_with_data(vec![vec![1.0; 2]; 2], points).is_err());
                }
                let mut weights = vec![vec![1.0; 2]; 2];
                weights[i][j] = invalid;
                assert!(surface_with_data(weights, vec![vec![vec![0.0]; 2]; 2]).is_err());
            }
        }
    }
    for weights in [vec![vec![0.0; 2]; 2], vec![vec![1.0, 1.0], vec![1.0, -1.0]]] {
        assert!(surface_with_data(weights, vec![vec![vec![0.0]; 2]; 2]).is_err());
    }
}

#[test]
fn surface_constructor_validates_both_knot_vectors_and_degrees() {
    let valid = vec![0.0, 0.0, 1.0, 1.0];
    for invalid_knots in [
        vec![],
        vec![0.0, 1.0],
        vec![1.0; 4],
        vec![0.0, 0.7, 0.2, 1.0],
        vec![0.0, 0.0, 0.0, 1.0],
        vec![f64::NAN, 0.0, 1.0, 1.0],
        vec![f64::NEG_INFINITY, 0.0, 1.0, 1.0],
        vec![0.0, 0.0, 1.0, f64::INFINITY],
    ] {
        for (u, v) in [
            (invalid_knots.clone(), valid.clone()),
            (valid.clone(), invalid_knots),
        ] {
            assert!(NURBSSurface::new(
                1,
                1,
                u,
                v,
                vec![vec![1.0; 2]; 2],
                vec![vec![vec![0.0]; 2]; 2]
            )
            .is_err());
        }
    }
    for invalid_degree in [2, usize::MAX - 1, usize::MAX] {
        for (u, v) in [(invalid_degree, 1), (1, invalid_degree)] {
            assert!(NURBSSurface::new(
                u,
                v,
                valid.clone(),
                valid.clone(),
                vec![vec![1.0; 2]; 2],
                vec![vec![vec![0.0]; 2]; 2]
            )
            .is_err());
        }
    }
}

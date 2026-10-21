use bsfun::{
    bspline_basis,
    nurbs::{NURBSCurve, NURBSSurface},
    nurbs_curve_point, nurbs_surface_point, rational_bspline_basis_curve,
    rational_bspline_basis_surface,
};

fn close(actual: f64, expected: f64) {
    assert!(actual.is_finite(), "non-finite result: {actual}");
    assert!(
        (actual - expected).abs() <= 2e-13 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

fn point_close(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len());
    for (&actual, &expected) in actual.iter().zip(expected) {
        close(actual, expected);
    }
}

fn clamped_knots(degree: usize) -> Vec<f64> {
    vec![0.0; degree + 1]
        .into_iter()
        .chain(vec![1.0; degree + 1])
        .collect()
}

#[test]
fn cubic_bases_match_bernstein_polynomials_including_endpoints() {
    let knots = clamped_knots(3);
    for sample in 0..=100 {
        let t = f64::from(sample) / 100.0;
        let expected = [
            (1.0 - t).powi(3),
            3.0 * t * (1.0 - t).powi(2),
            3.0 * t * t * (1.0 - t),
            t.powi(3),
        ];
        for (i, value) in expected.into_iter().enumerate() {
            close(bspline_basis(i, 3, &knots, t), value);
        }
    }
}

#[test]
fn nonuniform_bases_are_nonnegative_local_and_partition_unity() {
    for degree in 0..=6 {
        let knots: Vec<_> = vec![-2.0; degree + 1]
            .into_iter()
            .chain([-1.7, -0.25, 0.5, 2.75])
            .chain(vec![3.0; degree + 1])
            .collect();
        let count = knots.len() - degree - 1;
        for sample in 0..=500 {
            let t = -2.0 + f64::from(sample) / 100.0;
            let values: Vec<_> = (0..count)
                .map(|i| bspline_basis(i, degree, &knots, t))
                .collect();
            close(values.iter().sum(), 1.0);
            for (i, &value) in values.iter().enumerate() {
                assert!((0.0..=1.0 + 1e-14).contains(&value));
                if t < knots[i] || t > knots[i + degree + 1] {
                    assert_eq!(value, 0.0);
                }
            }
        }
    }
}

#[test]
fn repeated_interior_knots_use_the_right_hand_value() {
    let knots = [0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0];
    for i in 0..6 {
        close(
            bspline_basis(i, 2, &knots, 0.5),
            if i == 3 { 1.0 } else { 0.0 },
        );
    }
    for t in [0.0, 0.25, 0.5, 0.75, 1.0] {
        close((0..6).map(|i| bspline_basis(i, 2, &knots, t)).sum(), 1.0);
    }
}

#[test]
fn unclamped_right_endpoint_uses_left_limit_even_when_repeated() {
    let knots = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
    point_close(
        &(0..4)
            .map(|i| bspline_basis(i, 2, &knots, 4.0))
            .collect::<Vec<_>>(),
        &[0.0, 0.0, 0.5, 0.5],
    );
    let knots = [0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 5.0];
    point_close(
        &(0..4)
            .map(|i| bspline_basis(i, 2, &knots, 4.0))
            .collect::<Vec<_>>(),
        &[0.0, 0.0, 0.0, 1.0],
    );
}

#[test]
fn degree_zero_curves_and_surfaces_are_piecewise_constant() {
    let points = vec![vec![2.0], vec![7.0]];
    let curve = NURBSCurve::new(0, vec![0.0, 0.5, 1.0], vec![1.0, 3.0], points).unwrap();
    for (t, expected) in [(0.0, 2.0), (0.49, 2.0), (0.5, 7.0), (1.0, 7.0)] {
        point_close(&curve.eval(t).unwrap(), &[expected]);
    }
    let surface = NURBSSurface::new(
        0,
        0,
        vec![0.0, 1.0],
        vec![0.0, 1.0],
        vec![vec![2.0]],
        vec![vec![vec![3.0, 4.0]]],
    )
    .unwrap();
    point_close(&surface.eval(1.0, 1.0).unwrap(), &[3.0, 4.0]);
}

#[test]
fn rational_quarter_circle_has_known_midpoint_and_unit_radius() {
    let knots = clamped_knots(2);
    let weights = vec![1.0, std::f64::consts::FRAC_1_SQRT_2, 1.0];
    let points = vec![vec![1.0, 0.0], vec![1.0, 1.0], vec![0.0, 1.0]];
    let curve = NURBSCurve::new(2, knots.clone(), weights.clone(), points.clone()).unwrap();
    point_close(
        &curve.eval(0.5).unwrap(),
        &[std::f64::consts::FRAC_1_SQRT_2; 2],
    );
    for sample in 0..=100 {
        let t = f64::from(sample) / 100.0;
        let point = curve.eval(t).unwrap();
        close(point[0].hypot(point[1]), 1.0);
        point_close(&point, &nurbs_curve_point(t, &points, &weights, &knots, 2));
        close(
            (0..3)
                .map(|i| rational_bspline_basis_curve(i, 2, &knots, &weights, t))
                .sum(),
            1.0,
        );
    }
    point_close(&curve.eval(0.0).unwrap(), &points[0]);
    point_close(&curve.eval(1.0).unwrap(), &points[2]);
}

fn surface_points() -> Vec<Vec<Vec<f64>>> {
    vec![
        vec![vec![0.0, 0.0, 0.0], vec![0.0, 1.0, 0.0]],
        vec![vec![1.0, 0.0, 0.0], vec![1.0, 1.0, 1.0]],
    ]
}

#[test]
fn unequal_surface_weights_are_applied_once_and_bases_partition_unity() {
    let knots = clamped_knots(1);
    let weights = vec![vec![1.0, 2.0], vec![3.0, 4.0]];
    let points = surface_points();
    let surface = NURBSSurface::new(
        1,
        1,
        knots.clone(),
        knots.clone(),
        weights.clone(),
        points.clone(),
    )
    .unwrap();
    // At the centre all polynomial bases are 1/4, so rational bases are w/10.
    point_close(&surface.eval(0.5, 0.5).unwrap(), &[0.7, 0.6, 0.4]);
    for a in 0..=10 {
        for b in 0..=10 {
            let u = f64::from(a) / 10.0;
            let v = f64::from(b) / 10.0;
            let polynomial = [
                [(1.0 - u) * (1.0 - v), (1.0 - u) * v],
                [u * (1.0 - v), u * v],
            ];
            let denominator = polynomial[0][0]
                + 2.0 * polynomial[0][1]
                + 3.0 * polynomial[1][0]
                + 4.0 * polynomial[1][1];
            let mut sum = 0.0;
            for i in 0..2 {
                for j in 0..2 {
                    let basis =
                        rational_bspline_basis_surface(i, j, 1, 1, &knots, &knots, &weights, u, v);
                    close(basis, polynomial[i][j] * weights[i][j] / denominator);
                    sum += basis;
                }
            }
            close(sum, 1.0);
            let point = surface.eval(u, v).unwrap();
            point_close(
                &point,
                &[
                    (3.0 * polynomial[1][0] + 4.0 * polynomial[1][1]) / denominator,
                    (2.0 * polynomial[0][1] + 4.0 * polynomial[1][1]) / denominator,
                    4.0 * polynomial[1][1] / denominator,
                ],
            );
            point_close(
                &point,
                &nurbs_surface_point(&points, &weights, &knots, &knots, 1, 1, u, v),
            );
        }
    }
    for (i, row) in points.iter().enumerate() {
        for (j, point) in row.iter().enumerate() {
            point_close(&surface.eval(i as f64, j as f64).unwrap(), point);
        }
    }
}

#[test]
fn unit_weight_surface_is_bilinear_in_arbitrary_dimension() {
    let knots = clamped_knots(1);
    let points = vec![
        vec![vec![0.0, 0.0, 0.0, 5.0], vec![1.0, 0.0, 0.0, 5.0]],
        vec![vec![0.0, 1.0, 0.0, 5.0], vec![1.0, 1.0, 0.0, 5.0]],
    ];
    let surface =
        NURBSSurface::new(1, 1, knots.clone(), knots, vec![vec![1.0; 2]; 2], points).unwrap();
    point_close(&surface.eval(0.5, 0.75).unwrap(), &[0.75, 0.5, 0.0, 5.0]);
}

#[test]
fn common_weight_scaling_handles_tiny_huge_and_subnormal_weights() {
    let knots = clamped_knots(1);
    let points = vec![vec![2.0, -3.0], vec![6.0, 5.0]];
    for scale in [f64::from_bits(1), 1e-300, 1.0, 1e300, f64::MAX / 4.0] {
        let weights = vec![scale, 3.0 * scale];
        point_close(
            &nurbs_curve_point(0.5, &points, &weights, &knots, 1),
            &[5.0, 3.0],
        );
        let curve = NURBSCurve::new(1, knots.clone(), weights, points.clone()).unwrap();
        point_close(&curve.eval(0.5).unwrap(), &[5.0, 3.0]);
        let weights = vec![vec![scale, 2.0 * scale], vec![3.0 * scale, 4.0 * scale]];
        let surface = NURBSSurface::new(
            1,
            1,
            knots.clone(),
            knots.clone(),
            weights,
            surface_points(),
        )
        .unwrap();
        point_close(&surface.eval(0.5, 0.5).unwrap(), &[0.7, 0.6, 0.4]);
    }
}

#[test]
fn inactive_large_weights_do_not_underflow_supporting_weights() {
    let knots = vec![0.0, 0.0, 1.0, 2.0, 2.0];
    let curve = NURBSCurve::new(
        1,
        knots,
        vec![1e-300, 1e-300, f64::MAX],
        vec![vec![2.0], vec![4.0], vec![9.0]],
    )
    .unwrap();
    point_close(&curve.eval(0.5).unwrap(), &[3.0]);
    let surface = NURBSSurface::new(
        1,
        1,
        clamped_knots(1),
        clamped_knots(1),
        vec![vec![1e-300, f64::MAX], vec![1e-300, f64::MAX]],
        surface_points(),
    )
    .unwrap();
    point_close(&surface.eval(0.5, 0.0).unwrap(), &[0.5, 0.0, 0.0]);
}

#[test]
fn extreme_finite_coordinates_stay_finite_and_preserve_cancellation() {
    let knots = clamped_knots(1);
    for second in [f64::MAX, -f64::MAX] {
        let points = vec![vec![f64::MAX, 1e-300], vec![second, 1e-300]];
        let curve = NURBSCurve::new(1, knots.clone(), vec![1.0; 2], points).unwrap();
        let result = curve.eval(0.5).unwrap();
        close(result[0], if second > 0.0 { f64::MAX } else { 0.0 });
        assert_eq!(result[1], 1e-300);
    }
}

#[test]
fn wide_and_narrow_knot_intervals_avoid_overflow_and_epsilon_cutoffs() {
    for extent in [f64::from_bits(1), 1e-300, 1.0, 1e300] {
        let knots = [-extent, -extent, extent, extent];
        close(bspline_basis(0, 1, &knots, 0.0), 0.5);
        close(bspline_basis(1, 1, &knots, 0.0), 0.5);
    }
    let knots = [-f64::MAX, -f64::MAX, f64::MAX, f64::MAX];
    close(bspline_basis(0, 1, &knots, 0.0), 0.5);
    close(bspline_basis(1, 1, &knots, 0.0), 0.5);
}

#[test]
fn zero_denominators_are_errors_for_structs_and_zero_for_free_functions() {
    let knots = clamped_knots(1);
    let points = vec![vec![2.0], vec![4.0]];
    let weights = vec![0.0, 1.0];
    let curve = NURBSCurve::new(1, knots.clone(), weights.clone(), points.clone()).unwrap();
    assert!(curve.eval(0.0).is_err());
    point_close(&curve.eval(0.5).unwrap(), &[4.0]);
    point_close(
        &nurbs_curve_point(0.0, &points, &weights, &knots, 1),
        &[0.0],
    );
    assert_eq!(
        rational_bspline_basis_curve(0, 1, &knots, &weights, 0.0),
        0.0
    );
    let weights = vec![vec![0.0, 0.0], vec![1.0, 1.0]];
    let surface = NURBSSurface::new(
        1,
        1,
        knots.clone(),
        knots.clone(),
        weights.clone(),
        surface_points(),
    )
    .unwrap();
    assert!(surface.eval(0.0, 0.5).is_err());
    assert_eq!(
        rational_bspline_basis_surface(1, 1, 1, 1, &knots, &knots, &weights, 0.0, 0.5),
        0.0
    );
    point_close(
        &nurbs_surface_point(&surface_points(), &weights, &knots, &knots, 1, 1, 0.0, 0.5),
        &[0.0; 3],
    );
}

#[test]
fn struct_evaluation_rejects_nonfinite_and_outside_actual_domain() {
    let knots = vec![0.0, 1.0, 2.0, 3.0];
    let curve =
        NURBSCurve::new(1, knots.clone(), vec![1.0; 2], vec![vec![2.0], vec![4.0]]).unwrap();
    let surface = NURBSSurface::new(
        1,
        1,
        knots.clone(),
        knots.clone(),
        vec![vec![1.0; 2]; 2],
        surface_points(),
    )
    .unwrap();
    for t in [
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
        0.0,
        0.99,
        2.01,
        3.0,
    ] {
        assert!(curve.eval(t).is_err());
        assert!(surface.eval(t, 1.5).is_err());
        assert!(surface.eval(1.5, t).is_err());
        if !t.is_finite() {
            assert_eq!(bspline_basis(0, 1, &knots, t), 0.0);
        }
    }
    point_close(&curve.eval(1.0).unwrap(), &[2.0]);
    point_close(&curve.eval(2.0).unwrap(), &[4.0]);
}

#[test]
fn high_degree_point_evaluation_uses_the_private_iterative_helper() {
    let degree = 64;
    let knots = clamped_knots(degree);
    let weights = vec![1.0; degree + 1];
    let curve =
        NURBSCurve::new(degree, knots.clone(), weights, vec![vec![3.0]; degree + 1]).unwrap();
    point_close(&curve.eval(0.5).unwrap(), &[3.0]);
    let surface = NURBSSurface::new(
        degree,
        1,
        knots,
        clamped_knots(1),
        vec![vec![1.0; 2]; degree + 1],
        vec![vec![vec![3.0]; 2]; degree + 1],
    )
    .unwrap();
    point_close(&surface.eval(0.5, 0.5).unwrap(), &[3.0]);
}

#[test]
fn disparate_weights_are_scaled_by_their_actual_contributions() {
    let knots = clamped_knots(1);
    let weights = vec![1e-300, 1e300];
    let points = vec![vec![1e300], vec![0.0]];
    let curve = NURBSCurve::new(1, knots.clone(), weights.clone(), points).unwrap();
    // Weighted bases are approximately [1e-300, 1], so the first point still
    // contributes 1 to the result. Dividing weights by 1e300 first loses it.
    close(curve.eval(1e-300).unwrap()[0], 1.0);
    let basis = rational_bspline_basis_curve(0, 1, &knots, &weights, 1e-300);
    assert!((basis / 1e-300 - 1.0).abs() < 1e-14);
}

#[test]
fn rectangular_surface_with_different_degrees_is_an_exact_cylinder_patch() {
    let circle_points = [[1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let circle_weights = [1.0, std::f64::consts::FRAC_1_SQRT_2, 1.0];
    let weights = circle_weights.iter().map(|w| vec![*w, 2.0 * w]).collect();
    let points = circle_points
        .iter()
        .map(|p| vec![vec![p[0], p[1], 0.0, 7.0], vec![p[0], p[1], 1.0, 7.0]])
        .collect();
    let surface =
        NURBSSurface::new(2, 1, clamped_knots(2), clamped_knots(1), weights, points).unwrap();
    for a in 0..=20 {
        for b in 0..=20 {
            let u = f64::from(a) / 20.0;
            let v = f64::from(b) / 20.0;
            let result = surface.eval(u, v).unwrap();
            close(result[0].hypot(result[1]), 1.0);
            close(result[2], 2.0 * v / (1.0 + v));
            close(result[3], 7.0);
        }
    }
    point_close(
        &surface.eval(0.5, 0.5).unwrap(),
        &[
            std::f64::consts::FRAC_1_SQRT_2,
            std::f64::consts::FRAC_1_SQRT_2,
            2.0 / 3.0,
            7.0,
        ],
    );
}

#[test]
fn small_coordinate_contributions_survive_large_cancelling_terms() {
    let curve = NURBSCurve::new(
        2,
        clamped_knots(2),
        vec![1.0, 0.5, 1.0],
        vec![vec![1e16], vec![1.0], vec![-1e16]],
    )
    .unwrap();
    // All three rational bases are 1/3 at the midpoint.
    close(curve.eval(0.5).unwrap()[0], 1.0 / 3.0);
    let surface = NURBSSurface::new(
        1,
        1,
        clamped_knots(1),
        clamped_knots(1),
        vec![vec![1.0; 2]; 2],
        vec![vec![vec![1e16], vec![1.0]], vec![vec![-1e16], vec![1.0]]],
    )
    .unwrap();
    close(surface.eval(0.5, 0.5).unwrap()[0], 0.5);
}

#[test]
fn raw_polynomial_bases_keep_full_support_beyond_the_struct_domain() {
    let knots = [0.0, 1.0, 2.0, 3.0];
    // Cox-de Boor defines these tails even though a curve's standard domain
    // for this knot vector is [1, 2]. Its polynomial bases sum to less than one
    // on those tails, so checked curve evaluation intentionally rejects them.
    close(bspline_basis(0, 1, &knots, 0.5), 0.5);
    close(bspline_basis(1, 1, &knots, 2.5), 0.5);
    for t in [-1.0, 4.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(bspline_basis(0, 1, &knots, t), 0.0);
        assert_eq!(bspline_basis(1, 1, &knots, t), 0.0);
    }
    let points = vec![vec![2.0], vec![4.0]];
    point_close(
        &nurbs_curve_point(0.5, &points, &[1.0, 1.0], &knots, 1),
        &[2.0],
    );
    point_close(
        &nurbs_curve_point(2.5, &points, &[1.0, 1.0], &knots, 1),
        &[4.0],
    );
}

#[test]
fn recursive_and_iterative_curves_agree_with_repeated_and_unclamped_knots() {
    let mut cases = vec![
        (2, vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]),
        (2, vec![0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 5.0]),
    ];
    for degree in 0..=5 {
        let interior = if degree == 0 {
            vec![-1.0, 0.25, 0.5]
        } else {
            vec![-1.0, 0.25, 0.25, 0.5]
        };
        let knots = vec![-2.0; degree + 1]
            .into_iter()
            .chain(interior)
            .chain(vec![2.0; degree + 1])
            .collect();
        cases.push((degree, knots));
    }
    for (degree, knots) in cases {
        let count = knots.len() - degree - 1;
        let weights: Vec<_> = (0..count)
            .map(|i| {
                if i % 4 == 2 {
                    0.0
                } else {
                    (i + 1) as f64 / 2.0
                }
            })
            .collect();
        // With unit-vector control points, the evaluated point contains every
        // rational basis coefficient. This compares the public basis functions
        // as well as the two point evaluation paths.
        let points: Vec<Vec<f64>> = (0..count)
            .map(|i| {
                (0..count)
                    .map(|axis| if i == axis { 1.0 } else { 0.0 })
                    .collect()
            })
            .collect();
        let curve =
            NURBSCurve::new(degree, knots.clone(), weights.clone(), points.clone()).unwrap();
        let start = knots[degree];
        let end = knots[count];
        let parameters = knots
            .iter()
            .copied()
            .chain((0..=32).map(|sample| start + (end - start) * f64::from(sample) / 32.0))
            .filter(|t| *t >= start && *t <= end);
        for t in parameters {
            let recursive = nurbs_curve_point(t, &points, &weights, &knots, degree);
            let rational: Vec<_> = (0..count)
                .map(|i| rational_bspline_basis_curve(i, degree, &knots, &weights, t))
                .collect();
            point_close(&recursive, &rational);
            match curve.eval(t) {
                Ok(iterative) => {
                    point_close(&recursive, &iterative);
                    close(rational.iter().sum(), 1.0);
                }
                Err(_) => assert!(recursive.iter().all(|&value| value == 0.0)),
            }
        }
    }
}

#[test]
fn recursive_and_iterative_surfaces_agree_for_mixed_degrees_and_rectangular_nets() {
    let mut cases = vec![(
        2,
        1,
        vec![0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 5.0],
        vec![-1.0, 0.0, 0.5, 1.5, 2.0],
    )];
    for (degree_u, degree_v) in [(0, 0), (0, 2), (2, 0), (1, 2), (2, 3)] {
        let interior_u = if degree_u == 0 {
            vec![0.25, 0.6]
        } else {
            vec![0.25, 0.25, 0.6]
        };
        let knots_u = vec![0.0; degree_u + 1]
            .into_iter()
            .chain(interior_u)
            .chain(vec![1.0; degree_u + 1])
            .collect();
        let knots_v = vec![0.0; degree_v + 1]
            .into_iter()
            .chain([0.2, 0.7, 0.9])
            .chain(vec![1.0; degree_v + 1])
            .collect();
        cases.push((degree_u, degree_v, knots_u, knots_v));
    }
    for (degree_u, degree_v, knots_u, knots_v) in cases {
        let rows = knots_u.len() - degree_u - 1;
        let columns = knots_v.len() - degree_v - 1;
        let weights: Vec<Vec<f64>> = (0..rows)
            .map(|i| {
                (0..columns)
                    .map(|j| {
                        if (i + j) % 4 == 2 {
                            0.0
                        } else {
                            (1 + i * columns + j) as f64
                        }
                    })
                    .collect()
            })
            .collect();
        let points: Vec<Vec<Vec<f64>>> = (0..rows)
            .map(|i| {
                (0..columns)
                    .map(|j| vec![i as f64, -(j as f64), (i * j) as f64, 7.0])
                    .collect()
            })
            .collect();
        let surface = NURBSSurface::new(
            degree_u,
            degree_v,
            knots_u.clone(),
            knots_v.clone(),
            weights.clone(),
            points.clone(),
        )
        .unwrap();
        let start_u = knots_u[degree_u];
        let end_u = knots_u[rows];
        let start_v = knots_v[degree_v];
        let end_v = knots_v[columns];
        let parameters_u: Vec<_> = knots_u
            .iter()
            .copied()
            .chain((0..=8).map(|sample| start_u + (end_u - start_u) * f64::from(sample) / 8.0))
            .filter(|t| *t >= start_u && *t <= end_u)
            .collect();
        let parameters_v: Vec<_> = knots_v
            .iter()
            .copied()
            .chain((0..=8).map(|sample| start_v + (end_v - start_v) * f64::from(sample) / 8.0))
            .filter(|t| *t >= start_v && *t <= end_v)
            .collect();
        for &u in &parameters_u {
            for &v in &parameters_v {
                let recursive = nurbs_surface_point(
                    &points, &weights, &knots_u, &knots_v, degree_u, degree_v, u, v,
                );
                let mut basis_sum = 0.0;
                let mut from_rational_bases = vec![0.0; 4];
                for (i, row) in points.iter().enumerate() {
                    for (j, point) in row.iter().enumerate() {
                        let basis = rational_bspline_basis_surface(
                            i, j, degree_u, degree_v, &knots_u, &knots_v, &weights, u, v,
                        );
                        basis_sum += basis;
                        for (coordinate, value) in from_rational_bases.iter_mut().zip(point) {
                            *coordinate += basis * value;
                        }
                    }
                }
                point_close(&recursive, &from_rational_bases);
                match surface.eval(u, v) {
                    Ok(iterative) => {
                        point_close(&recursive, &iterative);
                        close(basis_sum, 1.0);
                    }
                    Err(_) => {
                        assert_eq!(basis_sum, 0.0);
                        assert!(recursive.iter().all(|&value| value == 0.0));
                    }
                }
            }
        }
    }
}

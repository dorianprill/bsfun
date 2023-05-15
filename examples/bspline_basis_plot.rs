use bsfun::bspline_basis;
use plotters::{prelude::*, style::full_palette::ORANGE};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let colors = [RED, GREEN, BLUE, ORANGE, CYAN, MAGENTA, BLACK];

    for degree in 0..4 {
        // Create a knot vector with degree + 1 repeated values at the start and end
        let knots = vec![0.0; degree + 1]
            .into_iter()
            .chain((1..=degree).map(|i| i as f64))
            .chain(vec![degree as f64 + 1.0; degree + 1])
            .collect::<Vec<_>>();
        // TODO check if this version makes more sense for degree 1
        // let num_internal_knots = degree + 1;
        // let knots = vec![0.0; degree+1]
        // .into_iter()
        // .chain((1..=num_internal_knots).map(|i| i as f64))
        // .chain(vec![degree as f64 + 1.0; degree+1])
        // .collect::<Vec<_>>();

        let filename = format!("basis_degree_{}.png", degree);
        let root = BitMapBackend::new(&filename, (400, 300)).into_drawing_area();
        root.fill(&WHITE)?;

        let mut chart = ChartBuilder::on(&root)
            .caption(
                format!("B-Spline Basis Degree {}", degree),
                ("sans-serif", 20).into_font(),
            )
            .margin(5)
            .x_label_area_size(30)
            .y_label_area_size(30)
            .build_cartesian_2d(0f64..(degree + 2) as f64, 0f64..1.05f64)?;

        chart.configure_mesh().draw()?;

        for i in 0..=knots.len() - degree - 2 {
            // Change here
            chart.draw_series(LineSeries::new(
                (0..=100).map(|x| {
                    let t = x as f64 * (degree + 2) as f64 / 100.0;
                    (t, bspline_basis(i, degree, &knots, t))
                }),
                colors[i % colors.len()], // Use the index to select the color
            ))?;
        }
    }

    Ok(())
}

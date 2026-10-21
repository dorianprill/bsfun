use bsfun::nurbs::NURBSCurve;
use kuva::plot::{LabelStyle, LegendEntry, LegendShape, LinePlot, ScatterPlot};
use kuva::render::{layout::Layout, plots::Plot, render::collect_legend_entries};
use std::f64::consts::FRAC_1_SQRT_2;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let degree = 2;
    let knots = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
    let weights = vec![1.0, FRAC_1_SQRT_2, 1.0];
    let control_points = vec![vec![1.0, 0.0], vec![1.0, 1.0], vec![0.0, 1.0]];
    let controls: Vec<_> = control_points.iter().map(|p| (p[0], p[1])).collect();
    let knot_label = knots
        .iter()
        .map(f64::to_string)
        .collect::<Vec<_>>()
        .join(", ");
    let curve = NURBSCurve::new(degree, knots, weights, control_points)?;
    let samples = (0..=200)
        .map(|sample| {
            let t = f64::from(sample) / 200.0;
            curve.eval(t).map(|point| (point[0], point[1]))
        })
        .collect::<Result<Vec<_>, _>>()?;

    let plots: Vec<Plot> = vec![
        LinePlot::new()
            .with_data(samples)
            .with_color("#0072B2")
            .with_stroke_width(3.0)
            .with_legend("NURBS curve")
            .into(),
        LinePlot::new()
            .with_data(controls.iter().copied())
            .with_color("#D55E00")
            .with_dashed()
            .with_legend("Control polygon")
            .into(),
        ScatterPlot::new()
            .with_data(controls)
            .with_color("#333333")
            .with_size(5.0)
            .with_labels(["P0", "P1", "P2"])
            .with_label_style(LabelStyle::Arrow {
                offset_x: 12.0,
                offset_y: -14.0,
            })
            .with_legend("Control points")
            .into(),
    ];
    let parameters = [
        format!("Degree = {degree}"),
        format!("Knots = [{knot_label}]"),
        "Weights = [1, 1/√2, 1]".into(),
        "Parameter t ∈ [0, 1]".into(),
    ]
    .into_iter()
    .map(|label| LegendEntry {
        label,
        // Parameter rows have text only, without a series swatch.
        color: "none".into(),
        shape: LegendShape::Line,
        dasharray: None,
    })
    .collect();
    let layout = Layout::auto_from_plots(&plots)
        .with_title("Exact NURBS quarter-circle")
        .with_x_label("x")
        .with_y_label("y")
        .with_x_axis_min(-0.12)
        .with_x_axis_max(1.15)
        .with_y_axis_min(-0.12)
        .with_y_axis_max(1.15)
        .with_equal_aspect()
        .with_width(900.0)
        .with_height(560.0)
        .with_legend_group("Geometry", collect_legend_entries(&plots))
        .with_legend_group("NURBS parameters", parameters);
    std::fs::write(
        "nurbs_quarter_circle.svg",
        kuva::render_to_svg(plots, layout),
    )?;
    println!("Wrote nurbs_quarter_circle.svg");
    Ok(())
}

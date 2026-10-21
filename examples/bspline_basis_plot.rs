use bsfun::bspline_basis;
use kuva::backend::svg::SvgBackend;
use kuva::plot::LinePlot;
use kuva::render::{figure::Figure, layout::Layout, plots::Plot};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let colors = [
        "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#333333",
    ];
    let mut panels = Vec::new();
    let mut layouts = Vec::new();

    for degree in 0..=3 {
        // Keep the same interior knots and domain in every panel so the effect
        // of degree is easy to compare. End knots have multiplicity degree + 1.
        let knots: Vec<_> = vec![0.0; degree + 1]
            .into_iter()
            .chain([1.0, 2.0, 3.0])
            .chain(vec![4.0; degree + 1])
            .collect();
        let plots: Vec<Plot> = (0..knots.len() - degree - 1)
            .map(|i| {
                let line = LinePlot::new()
                    .with_data((0..=400).map(|sample| {
                        let t = f64::from(sample) / 100.0;
                        (t, bspline_basis(i, degree, &knots, t))
                    }))
                    .with_color(colors[i % colors.len()])
                    .with_stroke_width(2.5)
                    .with_legend(format!("B{i}"));
                // Degree-zero bases jump at knots; don't draw sloping ramps.
                if degree == 0 { line.with_step() } else { line }.into()
            })
            .collect();
        layouts.push(
            Layout::new((0.0, 4.0), (0.0, 1.05))
                .with_title(format!("Degree {degree}"))
                .with_x_label("t")
                .with_y_label("Basis value"),
        );
        panels.push(plots);
    }

    let scene = Figure::new(2, 2)
        .with_plots(panels)
        .with_layouts(layouts)
        .with_title("B-spline bases · clamped knots, domain [0, 4]")
        .with_shared_legend_bottom()
        .with_cell_size(520.0, 340.0)
        .with_spacing(24.0)
        .render();
    std::fs::write("bspline_basis.svg", SvgBackend.render_scene(&scene))?;
    println!("Wrote bspline_basis.svg");
    Ok(())
}

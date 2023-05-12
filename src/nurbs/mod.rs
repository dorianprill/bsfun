// nurbs/mod.rs
pub mod nurbs_curve;
pub mod nurbs_surface;

// re-export to crate::nurbs::{...}
pub use nurbs_curve::NURBSCurve;
pub use nurbs_surface::NURBSSurface;

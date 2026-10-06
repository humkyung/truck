//! `IncludeCurve` on surfaces revolved from a line (an annulus or a cylinder made by `rsweep`).

use std::f64::consts::PI;
use truck_modeling::*;

/// Revolving a radial line gives a planar annulus sector; its own edges lie on it, a lifted edge does not.
/// (This used to recurse forever: the line entity was never converted.)
#[test]
fn include_on_revolved_line() {
    let a = builder::vertex(Point3::new(0.5, 0.0, 0.0));
    let b = builder::vertex(Point3::new(1.0, 0.0, 0.0));
    let line = builder::line(&a, &b);
    let shell: Shell =
        builder::rsweep(&line, Point3::origin(), Vector3::unit_z(), Rad(PI / 2.0), 1);
    let face = &shell[0];
    let surface = face.surface();
    assert!(matches!(surface, Surface::RevolutedCurve(_)));
    assert!(face.edge_iter().all(|edge| surface.include(&edge.curve())));
    let lifted = Curve::Line(Line(Point3::new(0.6, 0.0, 0.5), Point3::new(0.9, 0.0, 0.5)));
    assert!(!surface.include(&lifted));
}

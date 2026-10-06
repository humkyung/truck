//! Fillets on edges where at least one face is curved (plane–cylinder, plane–sphere, cylinder–cylinder):
//! exact rolling-ball surfaces, no panics, and `None` where no fillet exists.

use std::f64::consts::PI;
use std::panic::{catch_unwind, AssertUnwindSafe};
use truck_geometry::prelude::*;
use truck_modeling::{builder, Curve, Surface};
use truck_shapeops::fillet::{fillet_with_side, simple_fillet, FilletWithSide, SimpleFillet};
use truck_topology::shell::ShellCondition;

type Shell = truck_topology::Shell<Point3, Curve, Surface>;
type Solid = truck_topology::Solid<Point3, Curve, Surface>;

/// Index of the face whose surface passes through `p`.
fn face_at(shell: &Shell, p: Point3) -> usize {
    shell
        .face_iter()
        .position(|face| face.surface().search_parameter(p, None, 20).is_some())
        .unwrap_or_else(|| panic!("no face at {p:?}"))
}

/// The edge whose curve passes through `p` (inside its range).
fn edge_at(shell: &Shell, p: Point3) -> truck_topology::Edge<Point3, Curve> {
    shell
        .edge_iter()
        .find(|edge| {
            let curve = edge.curve();
            curve
                .search_nearest_parameter(p, None, 20)
                .is_some_and(|t| curve.subs(t).near(&p))
        })
        .unwrap_or_else(|| panic!("no edge at {p:?}"))
}

/// Result of a fillet attempt.
#[derive(Debug, PartialEq)]
enum Outcome {
    Panic,
    Refused,
    Filleted { closed: bool, consistent: bool },
}

/// Largest distance error of the fillet surface from the rolling-ball tube.
const TUBE_TOL: f64 = 1.0e-9;

/// The face order following the filleted edge gives a closed, consistent shell. The reversed order is the known trap
/// of `simple_fillet` (an open shell instead of `None`) — it must not panic.
fn check(outcome: Outcome, swap: bool) {
    match swap {
        false => assert_eq!(
            outcome,
            Outcome::Filleted {
                closed: true,
                consistent: true
            }
        ),
        true => assert_ne!(outcome, Outcome::Panic),
    }
}

/// Fillets the edge at `edge` between the faces at `f0`, `f1` (with the side faces at `sides`, if any) and checks the shell.
fn try_fillet(
    solid: &Solid,
    f0: Point3,
    f1: Point3,
    edge: Point3,
    sides: Option<(Point3, Point3)>,
    r: f64,
) -> (Outcome, Option<Shell>) {
    let mut shell = solid.boundaries()[0].clone();
    assert!(shell.is_geometric_consistent(), "input shell");
    let (i0, i1, e) = (
        face_at(&shell, f0),
        face_at(&shell, f1),
        edge_at(&shell, edge),
    );
    let sides = sides.map(|(s0, s1)| (face_at(&shell, s0), face_at(&shell, s1)));
    let run = || match sides {
        Some((s0, s1)) => fillet_with_side(
            &shell[i0],
            &shell[i1],
            e.id(),
            Some(&shell[s0]),
            Some(&shell[s1]),
            r,
            0.001,
        )
        .map(
            |FilletWithSide {
                 simple_fillet,
                 side0,
                 side1,
             }| (simple_fillet, side0, side1),
        ),
        None => simple_fillet(&shell[i0], &shell[i1], e.id(), r, 0.001).map(|s| (s, None, None)),
    };
    let result = match catch_unwind(AssertUnwindSafe(run)) {
        Err(_) => return (Outcome::Panic, None),
        Ok(None) => return (Outcome::Refused, None),
        Ok(Some(result)) => result,
    };
    let (
        SimpleFillet {
            face0,
            face1,
            fillet,
        },
        side0,
        side1,
    ) = result;
    shell[i0] = face0;
    shell[i1] = face1;
    if let (Some((s0, s1)), Some(n0), Some(n1)) = (sides, side0, side1) {
        shell[s0] = n0;
        shell[s1] = n1;
    }
    shell.push(fillet);
    let closed = shell.shell_condition() == ShellCondition::Closed;
    let outcome = Outcome::Filleted {
        closed,
        consistent: shell.is_geometric_consistent(),
    };
    (outcome, Some(shell))
}

/// Largest distance error of the fillet face (last face) from the tube of radius `r` around `centre_line`.
fn tube_error(shell: &Shell, r: f64, centre_distance: impl Fn(Point3) -> f64) -> f64 {
    let surface = shell.face_iter().last().unwrap().surface();
    let bound = |b: std::ops::Bound<f64>| match b {
        std::ops::Bound::Included(x) | std::ops::Bound::Excluded(x) => x,
        std::ops::Bound::Unbounded => panic!("unbounded fillet surface"),
    };
    let ((u0, u1), (v0, v1)) = surface.parameter_range();
    let ((u0, u1), (v0, v1)) = ((bound(u0), bound(u1)), (bound(v0), bound(v1)));
    let mut err: f64 = 0.0;
    for i in 0..=8 {
        for j in 0..=8 {
            let (u, v) = (
                u0 + (u1 - u0) * i as f64 / 8.0,
                v0 + (v1 - v0) * j as f64 / 8.0,
            );
            err = err.max((centre_distance(surface.subs(u, v)) - r).abs());
        }
    }
    err
}

/// Pie-slice prism: radius 1, quarter angle, height 1 (centre at the origin, slice between +x and +y).
fn pie_slice() -> Solid {
    let o = builder::vertex(Point3::origin());
    let a = builder::vertex(Point3::new(1.0, 0.0, 0.0));
    let b = builder::vertex(Point3::new(0.0, 1.0, 0.0));
    let transit = Point3::new(f64::cos(PI / 4.0), f64::sin(PI / 4.0), 0.0);
    let wire: truck_topology::Wire<Point3, Curve> = vec![
        builder::line(&o, &a),
        builder::circle_arc(&a, &b, transit),
        builder::line(&b, &o),
    ]
    .into();
    let face = builder::try_attach_plane(&[wire]).unwrap();
    builder::tsweep(&face, Vector3::unit_z())
}

/// Quarter of a spherical zone without the pole: the profile (bottom line, unit-sphere arc, vertical line at rho = 0.2)
/// revolved by 90° about z — faces: bottom plane, sphere, inner cylinder, two side planes.
fn sphere_quarter() -> Solid {
    let rho = 0.2;
    let (top_angle, z_top) = (f64::acos(rho), f64::sqrt(1.0 - rho * rho));
    let o = builder::vertex(Point3::new(rho, 0.0, 0.0));
    let a = builder::vertex(Point3::new(1.0, 0.0, 0.0));
    let top = builder::vertex(Point3::new(rho, 0.0, z_top));
    let mid = top_angle / 2.0;
    let wire: truck_topology::Wire<Point3, Curve> = vec![
        builder::line(&o, &a),
        builder::circle_arc(&a, &top, Point3::new(f64::cos(mid), 0.0, f64::sin(mid))),
        builder::line(&top, &o),
    ]
    .into();
    let face = builder::try_attach_plane(&[wire]).unwrap();
    builder::rsweep(&face, Point3::origin(), Vector3::unit_z(), Rad(PI / 2.0), 1)
}

/// Lens prism: intersection of two unit disks centred at (±a, 0), height 1 — the cusps at (0, ±a) are 90° corners.
fn lens() -> Solid {
    let a = f64::sqrt(0.5);
    let p = builder::vertex(Point3::new(0.0, -a, 0.0));
    let q = builder::vertex(Point3::new(0.0, a, 0.0));
    let wire: truck_topology::Wire<Point3, Curve> = vec![
        builder::circle_arc(&p, &q, Point3::new(1.0 - a, 0.0, 0.0)),
        builder::circle_arc(&q, &p, Point3::new(a - 1.0, 0.0, 0.0)),
    ]
    .into();
    let face = builder::try_attach_plane(&[wire]).unwrap();
    builder::tsweep(&face, Vector3::unit_z())
}

/// Full cylinder: radius 1, height 1, rims made of two half circles.
fn cylinder() -> Solid {
    let p = builder::vertex(Point3::new(1.0, 0.0, 0.0));
    let q = builder::vertex(Point3::new(-1.0, 0.0, 0.0));
    let wire: truck_topology::Wire<Point3, Curve> = vec![
        builder::circle_arc(&p, &q, Point3::new(0.0, 1.0, 0.0)),
        builder::circle_arc(&q, &p, Point3::new(0.0, -1.0, 0.0)),
    ]
    .into();
    let face = builder::try_attach_plane(&[wire]).unwrap();
    builder::tsweep(&face, Vector3::unit_z())
}

#[test]
fn plane_cylinder_rim() {
    let (c, s) = (f64::cos(PI / 4.0), f64::sin(PI / 4.0));
    let solid = pie_slice();
    for swap in [false, true] {
        let (top, side) = (Point3::new(0.3, 0.3, 1.0), Point3::new(c, s, 0.5));
        let (f0, f1) = if swap { (side, top) } else { (top, side) };
        let r = 0.1;
        let (outcome, shell) = try_fillet(
            &solid,
            f0,
            f1,
            Point3::new(c, s, 1.0),
            Some((Point3::new(0.5, 0.0, 0.5), Point3::new(0.0, 0.5, 0.5))),
            r,
        );
        check(outcome, swap);
        if let Some(shell) = shell {
            let err = tube_error(&shell, r, |p| {
                let rho = Vector2::new(p.x, p.y).magnitude();
                Vector2::new(rho - (1.0 - r), p.z - (1.0 - r)).magnitude()
            });
            assert!(err < TUBE_TOL, "tube error {err:e}");
        }
    }
}

#[test]
fn plane_sphere_rim() {
    let (c, s) = (f64::cos(PI / 4.0), f64::sin(PI / 4.0));
    let solid = sphere_quarter();
    for swap in [false, true] {
        let (bottom, ball) = (
            Point3::new(0.3, 0.3, 0.0),
            Point3::new(0.5, 0.5, f64::sqrt(0.5)),
        );
        let (f0, f1) = if swap { (ball, bottom) } else { (bottom, ball) };
        let r = 0.1;
        let (outcome, shell) = try_fillet(
            &solid,
            f0,
            f1,
            Point3::new(c, s, 0.0),
            Some((Point3::new(0.3, 0.0, 0.3), Point3::new(0.0, 0.3, 0.3))),
            r,
        );
        check(outcome, swap);
        if let Some(shell) = shell {
            // The ball centre runs on the circle rho = sqrt((1-r)^2 - r^2), z = r.
            let rho0 = f64::sqrt((1.0 - r) * (1.0 - r) - r * r);
            let err = tube_error(&shell, r, |p| {
                let rho = Vector2::new(p.x, p.y).magnitude();
                Vector2::new(rho - rho0, p.z - r).magnitude()
            });
            assert!(err < TUBE_TOL, "tube error {err:e}");
        }
    }
}

#[test]
fn cylinder_cylinder_edge() {
    let a = f64::sqrt(0.5);
    let solid = lens();
    for swap in [false, true] {
        let (right, left) = (
            Point3::new(1.0 - a, 0.0, 0.5),
            Point3::new(a - 1.0, 0.0, 0.5),
        );
        let (f0, f1) = if swap { (left, right) } else { (right, left) };
        let r = 0.1;
        let (outcome, shell) = try_fillet(
            &solid,
            f0,
            f1,
            Point3::new(0.0, a, 0.5),
            Some((Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 1.0))),
            r,
        );
        check(outcome, swap);
        if let Some(shell) = shell {
            // The ball centre is r inside both cylinders: distance 1 - r from both axes (±a, 0).
            let cy = f64::sqrt((1.0 - r) * (1.0 - r) - a * a);
            let err = tube_error(&shell, r, |p| Vector2::new(p.x, p.y - cy).magnitude());
            assert!(err < TUBE_TOL, "tube error {err:e}");
        }
    }
}

#[test]
fn full_cylinder_rim() {
    let solid = cylinder();
    for swap in [false, true] {
        let (top, side) = (Point3::new(0.0, 0.5, 1.0), Point3::new(0.0, 1.0, 0.5));
        let (f0, f1) = if swap { (side, top) } else { (top, side) };
        // The half circles meet tangentially: the contact curve never crosses the neighbouring edge.
        let (outcome, _) = try_fillet(&solid, f0, f1, Point3::new(0.0, 1.0, 1.0), None, 0.1);
        assert_eq!(outcome, Outcome::Refused);
    }
}

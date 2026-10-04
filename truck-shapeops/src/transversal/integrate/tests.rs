use truck_meshalgo::prelude::*;
use truck_modeling::*;
use truck_topology::shell::ShellCondition;

#[test]
fn punched_cube() {
    let v = builder::vertex(Point3::origin());
    let e = builder::tsweep(&v, Vector3::unit_x());
    let f = builder::tsweep(&e, Vector3::unit_y());
    let cube: Solid = builder::tsweep(&f, Vector3::unit_z());

    let v = builder::vertex(Point3::new(0.5, 0.25, -0.5));
    let w = builder::rsweep(
        &v,
        Point3::new(0.5, 0.5, 0.0),
        Vector3::unit_z(),
        Rad(7.0),
        3,
    );
    let f = builder::try_attach_plane(&[w]).unwrap();
    let mut cylinder = builder::tsweep(&f, Vector3::unit_z() * 2.0);
    cylinder.not();
    let and = crate::and(&cube, &cylinder, 0.05).unwrap();

    let poly = and.triangulation(0.01).to_polygon();
    let file = std::fs::File::create("punched-cube.obj").unwrap();
    obj::write(&poly, file).unwrap();
}

fn unit_cube() -> Solid {
    let v = builder::vertex(Point3::origin());
    let e = builder::tsweep(&v, Vector3::unit_x());
    let f = builder::tsweep(&e, Vector3::unit_y());
    builder::tsweep(&f, Vector3::unit_z())
}

/// A vertical cylinder with axis through `(cx, cy)`, radius `r`, spanning `z0 .. z0 + h`.
fn cylinder(cx: f64, cy: f64, r: f64, z0: f64, h: f64) -> Solid {
    let v = builder::vertex(Point3::new(cx + r, cy, z0));
    let w = builder::rsweep(&v, Point3::new(cx, cy, z0), Vector3::unit_z(), Rad(7.0), 3);
    let f = builder::try_attach_plane(&[w]).unwrap();
    builder::tsweep(&f, Vector3::unit_z() * h)
}

fn volume(solid: &Solid) -> f64 {
    let mut poly = solid.triangulation(0.01).to_polygon();
    poly.put_together_same_attrs(TOLERANCE)
        .remove_degenerate_faces()
        .remove_unused_attrs();
    assert_eq!(poly.shell_condition(), ShellCondition::Closed);
    let ps = poly.positions();
    poly.faces()
        .triangle_iter()
        .map(|t| {
            let [a, b, c] = [ps[t[0].pos], ps[t[1].pos], ps[t[2].pos]];
            a.to_vec().dot(b.to_vec().cross(c.to_vec())) / 6.0
        })
        .sum()
}

/// Difference with a tool that does not touch the target returns the target unchanged.
#[test]
fn difference_with_disjoint_tool() {
    let cube = unit_cube();
    let mut tool = cylinder(5.0, 5.0, 0.25, -0.5, 2.0);
    tool.not();
    let res = crate::and(&cube, &tool, 0.05).unwrap();
    assert_eq!(res.boundaries().len(), 1);
    assert_eq!(res.boundaries()[0].len(), 6);
    assert_near!(volume(&res), 1.0);
}

/// Difference with a tool lying inside an existing hole (touching no material) returns the target unchanged.
#[test]
fn difference_with_tool_inside_hole() {
    let mut hole = cylinder(0.5, 0.5, 0.25, -0.5, 2.0);
    hole.not();
    let punched = crate::and(&unit_cube(), &hole, 0.05).unwrap();
    let before = volume(&punched);
    let mut tool = cylinder(0.5, 0.5, 0.1, -0.5, 2.0);
    tool.not();
    let res = crate::and(&punched, &tool, 0.05).unwrap();
    assert_eq!(res.boundaries().len(), 1);
    assert_eq!(res.boundaries()[0].len(), punched.boundaries()[0].len());
    assert!((volume(&res) - before).abs() < 1.0e-6);
}

#[test]
fn adjacent_cubes_or() {
    let v = builder::vertex(Point3::origin());
    let e = builder::tsweep(&v, Vector3::unit_x());
    let f = builder::tsweep(&e, Vector3::unit_y());
    let cube: Solid = builder::tsweep(&f, Vector3::unit_z());

    let v = builder::vertex(Point3::new(0.5, 0.5, 1.0));
    let w = builder::tsweep(&v, Vector3::unit_x());
    let f = builder::tsweep(&w, Vector3::unit_y());
    let cube2: Solid = builder::tsweep(&f, Vector3::unit_z());

    let result = crate::or(&cube, &cube2, 0.05);
    assert!(
        result.is_some(),
        "Boolean OR of adjacent cubes should succeed"
    );
    let solid = result.unwrap();

    // The result should be a single solid (one shell)
    assert_eq!(solid.boundaries().len(), 1);

    // Triangulate to check volume and bounding box
    let poly = solid.triangulation(0.01).to_polygon();

    // The volume of two unit cubes sharing only a portion of their boundary should be 2.0
    assert_near!(poly.volume(), 2.0);

    // Check the center of gravity
    let homog = poly.center_of_gravity();
    assert_near!(homog.to_point(), Point3::new(0.75, 0.75, 1.0));

    // Check the bounding box
    let bbx = poly.bounding_box();
    assert_near!(bbx.min(), Point3::new(0.0, 0.0, 0.0));
    assert_near!(bbx.max(), Point3::new(1.5, 1.5, 2.0));

    // Optional: check face count.
    // Two cubes normally have 6 faces each.
    // After union, the internal shared part is removed.
    // We expect around 12 faces as analyzed.
    assert_eq!(solid.face_iter().count(), 12);
}

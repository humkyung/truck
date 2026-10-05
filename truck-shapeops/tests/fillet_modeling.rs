//! Fillets with the geometry types of `truck_modeling` (fork: `Curve::PCurve`, `Surface::ApproxFillet`).

use truck_geometry::prelude::*;

/// Fillets a cube edge with the geometry types of `truck_modeling`, which can now hold
/// fillet surfaces (`Surface::ApproxFillet`) and their contact curves (`Curve::PCurve`).
#[test]
fn fillet_cube_with_modeling_geometry() {
    use truck_modeling::{primitive, Curve as MCurve, Surface as MSurface};
    use truck_shapeops::fillet::{fillet_with_side, FilletWithSide, SimpleFillet};
    use truck_topology::shell::ShellCondition;
    type MShell = truck_topology::Shell<Point3, MCurve, MSurface>;

    let bbox = BoundingBox::from_iter([Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)]);
    let solid: truck_topology::Solid<Point3, MCurve, MSurface> = primitive::cuboid(bbox);
    let mut shell: MShell = solid.into_boundaries().pop().unwrap();
    let face_at = |shell: &MShell, p: Point3| {
        shell
            .face_iter()
            .position(|face| face.surface().search_parameter(p, None, 10).is_some())
            .unwrap()
    };
    let (face0, face1) = (
        face_at(&shell, (0.5, 0.0, 0.5).into()),
        face_at(&shell, (1.0, 0.5, 0.5).into()),
    );
    let (side0, side1) = (
        face_at(&shell, (0.5, 0.5, 0.0).into()),
        face_at(&shell, (0.5, 0.5, 1.0).into()),
    );
    let edge = shell
        .edge_iter()
        .find(|edge| {
            edge.curve()
                .search_parameter((1.0, 0.0, 0.5).into(), None, 10)
                .is_some()
        })
        .unwrap();
    let FilletWithSide {
        simple_fillet:
            SimpleFillet {
                fillet,
                face0: new0,
                face1: new1,
            },
        side0: new_side0,
        side1: new_side1,
    } = fillet_with_side(
        &shell[face0],
        &shell[face1],
        edge.id(),
        Some(&shell[side0]),
        Some(&shell[side1]),
        0.2,
        0.001,
    )
    .unwrap();
    assert!(matches!(fillet.surface(), MSurface::ApproxFillet(_)));
    shell[face0] = new0;
    shell[face1] = new1;
    shell[side0] = new_side0.unwrap();
    shell[side1] = new_side1.unwrap();
    shell.push(fillet);
    assert_eq!(shell.shell_condition(), ShellCondition::Closed);
    assert!(shell.is_geometric_consistent());
    // The fillet surface lies on the rolling ball: distance from the axis (1 - r, r, z) is r.
    let fillet = shell.face_iter().last().unwrap();
    let surface = fillet.surface();
    let axis = Point3::new(0.8, 0.2, 0.0);
    for p in [(0.1, 0.5), (0.5, 0.5), (0.9, 0.5)] {
        let q = ParametricSurface::subs(&surface, p.0, p.1);
        let d = Vector2::new(q.x - axis.x, q.y - axis.y).magnitude();
        assert!((d - 0.2).abs() < 1.0e-6, "{d}");
    }
    // Round trip through serde keeps the new variants.
    let json = serde_json::to_string(&shell.compress()).unwrap();
    let back: truck_topology::compress::CompressedShell<Point3, MCurve, MSurface> =
        serde_json::from_str(&json).unwrap();
    assert_eq!(back.faces.len(), shell.len());
}

#![allow(clippy::many_single_char_names)]

use super::faces_classification::FacesClassification;
use super::loops_store::*;
use rustc_hash::FxHashMap as HashMap;
use std::ops::Deref;
use truck_meshalgo::prelude::*;
use truck_topology::*;

fn create_parameter_boundary<P, C, S>(
    face: &Face<P, C, S>,
    wire: &Wire<P, C>,
    polys: &mut HashMap<EdgeID<C>, PolylineCurve<P>>,
    tol: f64,
) -> Option<PolylineCurve<Point2>>
where
    P: Copy,
    C: BoundedCurve<Point = P> + ParameterDivision1D<Point = P>,
    S: Clone + SearchParameter<D2, Point = P>,
{
    let surface = face.surface();
    let pt = wire.front_vertex().unwrap().point();
    let p: Point2 = surface.search_parameter(pt, None, 100)?.into();
    let vec = wire.edge_iter().try_fold(vec![p], |mut vec, edge| {
        let poly = polys.entry(edge.id()).or_insert_with(|| {
            let curve = edge.curve();
            let div = curve.parameter_division(curve.range_tuple(), tol).1;
            PolylineCurve(div)
        });
        let mut p = *vec.last().unwrap();
        let closure = |q: &P| -> Option<Point2> {
            p = surface.search_parameter(*q, Some(p.into()), 100)?.into();
            Some(p)
        };
        let add: Option<Vec<Point2>> = match edge.orientation() {
            true => poly.iter().skip(1).map(closure).collect(),
            false => poly.iter().rev().skip(1).map(closure).collect(),
        };
        vec.append(&mut add?);
        Some(vec)
    })?;
    Some(PolylineCurve(vec))
}

#[derive(Clone, Debug)]
struct WireChunk<'a, C> {
    poly: PolylineCurve<Point2>,
    wire: &'a BoundaryWire<Point3, C>,
}

// Coincident contours may have distinct edge IDs. Compare the entire sampled UV
// boundary, in opposite traversal order, up to floating-point roundoff. This is
// deliberately conservative: unequal sampling is not proof of cancellation.
fn same_parameter_boundary(a: &PolylineCurve<Point2>, b: &PolylineCurve<Point2>) -> bool {
    let scale = a
        .iter()
        .chain(b.iter())
        .flat_map(|p| [p.x.abs(), p.y.abs()])
        .fold(0.0, f64::max);
    let epsilon = 64.0 * f64::EPSILON * scale;
    let near = |p: Point2, q: Point2| (p.x - q.x).abs() <= epsilon && (p.y - q.y).abs() <= epsilon;
    let points = |poly: &PolylineCurve<Point2>| {
        let mut points = poly.iter().copied().collect::<Vec<_>>();
        if points.len() > 1 && near(points[0], *points.last().unwrap()) {
            points.pop();
        }
        points
    };
    let (a, b) = (points(a), points(b));
    a.len() >= 3
        && a.len() == b.len()
        && (0..b.len()).any(|start| {
            a.iter()
                .enumerate()
                .all(|(i, &p)| near(p, b[(start + b.len() - i) % b.len()]))
        })
}

type FaceWithShapesOpStatus<C, S> = (Face<Point3, C, S>, ShapesOpStatus);
fn divide_one_face<C, S>(
    face: &Face<Point3, C, S>,
    loops: &Loops<Point3, C>,
    tol: f64,
) -> Option<Vec<FaceWithShapesOpStatus<C, S>>>
where
    C: BoundedCurve<Point = Point3> + ParameterDivision1D<Point = Point3>,
    S: Clone + SearchParameter<D2, Point = Point3>,
{
    let (mut pre_faces, mut negative_wires) = (Vec::new(), Vec::new());
    let mut map = HashMap::default();
    loops.iter().try_for_each(|wire| {
        let poly = create_parameter_boundary(face, wire, &mut map, tol)?;
        let area = poly.area();
        if !area.is_finite() {
            return None;
        }
        // UV area and the model-space tessellation tolerance have different units.
        // Do not delete a small but nonzero loop: it may be a real hole or thin wall.
        if area == 0.0 {
            return Some(());
        }
        match area > 0.0 {
            true => pre_faces.push(vec![WireChunk { poly, wire }]),
            false => negative_wires.push(WireChunk { poly, wire }),
        }
        Some(())
    })?;
    negative_wires.into_iter().try_for_each(|chunk| {
        // Preserve containment order: a coincident positive loop may itself be a
        // separate face inside the actual owner. Do not cancel it globally.
        let i = pre_faces.iter().position(|face| face[0].poly.include(chunk.poly.front()))?;
        if same_parameter_boundary(&pre_faces[i][0].poly, &chunk.poly) {
            if pre_faces[i].len() != 1 {
                return None;
            }
            pre_faces.remove(i);
        } else {
            // Equal areas alone do not establish coincident boundaries.
            if -chunk.poly.area() >= pre_faces[i][0].poly.area() {
                return None;
            }
            pre_faces[i].push(chunk);
        }
        Some(())
    })?;
    let vec: Vec<_> = pre_faces
        .into_iter()
        .map(|pre_face| {
            let surface = face.surface();
            let op = pre_face
                .iter()
                .find(|chunk| chunk.wire.status() != ShapesOpStatus::Unknown);
            let status = match op {
                Some(chunk) => chunk.wire.status(),
                None => ShapesOpStatus::Unknown,
            };
            let wires: Vec<Wire<Point3, C>> = pre_face
                .into_iter()
                .map(|chunk| chunk.wire.deref().clone())
                .collect();
            let mut new_face = Face::try_new(wires, surface).ok()?;
            if !face.orientation() {
                new_face.invert();
            }
            Some((new_face, status))
        })
        .collect::<Option<_>>()?;
    Some(vec)
}

pub fn divide_faces<C, S>(
    shell: &Shell<Point3, C, S>,
    loops_store: &LoopsStore<Point3, C>,
    tol: f64,
) -> Option<FacesClassification<Point3, C, S>>
where
    C: BoundedCurve<Point = Point3> + ParameterDivision1D<Point = Point3>,
    S: Clone + SearchParameter<D2, Point = Point3>,
{
    let mut res = FacesClassification::<Point3, C, S>::default();
    shell
        .iter()
        .zip(loops_store)
        .try_for_each(|(face, loops)| {
            if loops
                .iter()
                .all(|wire| wire.status() == ShapesOpStatus::Unknown)
            {
                res.push(face.clone(), ShapesOpStatus::Unknown);
            } else {
                let vec = divide_one_face(face, loops, tol)?;
                vec.into_iter()
                    .for_each(|(face, status)| res.push(face, status));
            }
            Some(())
        })?;
    Some(res)
}

#[cfg(test)]
mod tests;

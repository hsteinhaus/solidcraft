//! The exact cut of a body by a sketch plane (`sketch.intersect`), face by face: a plane, a
//! cylinder along or across the plane's normal, a sphere, and a cone or torus about the normal
//! cut it in lines and circles, trimmed where the face's edges cross the plane (on its vertices
//! where the cut passes through them); edges lying in the plane are those edges themselves.
//! A spline face extruded straight along a direction in the plane is cut in its rulings through
//! the crossings of its edges. Any other face is cut on its mesh, as before.

use solidcore_record::brep::{Brep, Loop, Surface};
use solidcore_record::{Plane, Vec2, Vec3};
use solidcraft_geom::Mesh;
use solidcraft_geom::brep::eval::CurveEval;
use solidcraft_geom::prelude::*;
use solidcraft_sketch::LinkGeom;
use std::collections::BTreeSet;

use crate::eval::ModelBody;

/// Samples along an edge where its crossings with the plane are looked for.
const SAMPLES: usize = 64;
/// Sines below this are parallel (unit vectors).
const PARALLEL: f64 = 1e-12;

/// What the plane cuts from a face's surface, in plane coordinates.
enum Cut {
    /// Lines: a point and the unit direction.
    Lines(Vec<(Vec2, Vec2)>),
    /// Circles: centre and radius.
    Circles(Vec<(Vec2, f64)>),
    /// Lines of this unit direction through the crossings of the face's boundary: the rulings of
    /// a wall extruded along the plane.
    Rulings(Vec2),
    /// The face lies in the plane, or doesn't reach it.
    Nothing,
    /// Not an exact cut here.
    Unsupported,
}

/// A point where the plane cuts a face's boundary, and how near it must be to a cut curve.
#[derive(Clone, Copy)]
struct Crossing {
    at: Vec3,
    tol: f64,
}

/// `body` cut by `plane`, exact on analytic faces: `faces` only, when given. Other faces are cut
/// on `mesh` by `mesh_cut`. `None` for a body without a record.
pub(crate) fn section(
    body: &ModelBody,
    mesh: &Mesh,
    plane: &Plane,
    faces: Option<&[u32]>,
    mesh_cut: impl Fn(u32) -> Vec<LinkGeom>,
) -> Option<Vec<LinkGeom>> {
    if body.body.is_mesh() {
        return None;
    }
    let record = body.body.record();
    let size = body.body.size().max(1e-9);
    let tiny = size * 1e-9;
    let n = plane.normal().normalized()?;
    let shape = solidcraft_kernel::DistanceShape::new(&body.body, mesh, body.display_tol());
    let wanted = |f: u32| faces.is_none_or(|fs| fs.contains(&f));
    // Edges in the plane, and every edge's crossings with it.
    let mut flat: BTreeSet<usize> = BTreeSet::new();
    let mut crossings: Vec<Vec<Crossing>> = Vec::with_capacity(record.edges.len());
    for i in 0..record.edges.len() {
        let (lies, at) = edge_crossings(record, i, plane, tiny);
        if lies {
            flat.insert(i);
        }
        crossings.push(at.unwrap_or_default());
    }
    let mut out = Vec::new();
    // The edges in the plane, once each.
    let flat_wanted: Vec<usize> = flat
        .iter()
        .copied()
        .filter(|e| record.edges.get(*e).is_some_and(|x| if x.uses.is_empty() { faces.is_none() } else { x.uses.iter().any(|u| wanted(u.face)) }))
        .collect();
    for g in super::exact_edges(body, plane, &flat_wanted) {
        out.extend(g?);
    }
    for (fi, face) in record.faces.iter().enumerate() {
        let f = u32::try_from(fi).ok()?;
        if !wanted(f) {
            continue;
        }
        let cut = surface_cut(&face.surface, plane, n, tiny);
        if matches!(cut, Cut::Nothing) {
            continue;
        }
        if matches!(cut, Cut::Unsupported) {
            out.extend(mesh_cut(f));
            continue;
        }
        // The face's boundary crossings, its edges in the plane, and its seams (an edge the face
        // holds twice: the cut goes on across it).
        let mut uses: Vec<u32> = Vec::new();
        let mut points: Vec<Crossing> = Vec::new();
        for l in &face.loops {
            match l {
                Loop::Edges { edges } => uses.extend(edges.iter().map(|u| u.edge)),
                Loop::Vertex { vertex } => {
                    if let Some(v) = record.vertices.get(*vertex as usize) {
                        let at = Vec3::from(v.point);
                        let tol = v.tolerance.max(tiny);
                        if plane.height(at).abs() <= tol {
                            points.push(Crossing { at, tol });
                        }
                    }
                }
            }
        }
        let seams: BTreeSet<u32> = uses.iter().copied().filter(|e| uses.iter().filter(|x| *x == e).count() > 1).collect();
        let own_flat: Vec<usize> = uses.iter().map(|e| *e as usize).filter(|e| flat.contains(e)).collect();
        for e in uses.iter().filter(|e| !seams.contains(e)) {
            points.extend(crossings.get(*e as usize).into_iter().flatten().copied());
        }
        let on_flat = |p: Vec3| {
            own_flat.iter().any(|e| record.edges.get(*e).is_some_and(|x| distance_to_edge(x, p).is_some_and(|d| d <= x.tolerance.max(tiny))))
        };
        let inside = |p: Vec3| shape.on_face(f, p) && !on_flat(p);
        let cut = match cut {
            Cut::Rulings(d) => {
                let mut lines: Vec<(Vec2, Vec2)> = Vec::new();
                for c in &points {
                    let q = plane.to_local(c.at);
                    if !lines.iter().any(|(p, _)| (*p + d * (q - *p).dot(d)).dist(q) <= c.tol + tiny) {
                        lines.push((q, d));
                    }
                }
                // A wall the plane only touches has no crossing to tell where.
                if lines.is_empty() {
                    out.extend(mesh_cut(f));
                    continue;
                }
                Cut::Lines(lines)
            }
            cut => cut,
        };
        match cut {
            Cut::Lines(lines) => {
                for (p, d) in lines {
                    let mut at: Vec<(f64, Vec2)> = Vec::new();
                    for c in &points {
                        let q = plane.to_local(c.at);
                        let s = (q - p).dot(d);
                        if (p + d * s).dist(q) <= c.tol + tiny {
                            at.push((s, q));
                        }
                    }
                    at.sort_by(|a, b| a.0.total_cmp(&b.0));
                    at.dedup_by(|b, a| a.1.dist(b.1) <= tiny.max(1e-12));
                    for w in at.windows(2) {
                        let [(sa, a), (sb, b)] = [w[0], w[1]];
                        if sb - sa <= tiny {
                            continue;
                        }
                        if inside(plane.to_world(p + d * (0.5 * (sa + sb)))) {
                            out.push(LinkGeom::Line(a, b));
                        }
                    }
                }
            }
            Cut::Circles(circles) => {
                for (c, r) in circles {
                    let mut at: Vec<(f64, Vec2)> = Vec::new();
                    for x in &points {
                        let q = plane.to_local(x.at);
                        if (q.dist(c) - r).abs() <= x.tol + tiny {
                            at.push(((q - c).angle(), q));
                        }
                    }
                    at.sort_by(|a, b| a.0.total_cmp(&b.0));
                    at.dedup_by(|b, a| a.1.dist(b.1) <= tiny.max(1e-12));
                    if at.len() > 1 && at.first().zip(at.last()).is_some_and(|(a, b)| a.1.dist(b.1) <= tiny.max(1e-12)) {
                        at.pop();
                    }
                    let point = |t: f64| c + Vec2::from_angle(t) * r;
                    if at.is_empty() {
                        if inside(plane.to_world(point(0.0))) {
                            out.push(LinkGeom::Circle(c, r));
                        }
                        continue;
                    }
                    for i in 0..at.len() {
                        let (ta, a) = at[i];
                        let (tb, b) = at[(i + 1) % at.len()];
                        let sweep = (tb - ta).rem_euclid(std::f64::consts::TAU);
                        let sweep = if at.len() == 1 { std::f64::consts::TAU } else { sweep };
                        if sweep * r <= tiny {
                            continue;
                        }
                        if inside(plane.to_world(point(ta + 0.5 * sweep))) {
                            out.push(LinkGeom::Arc { c, a, b });
                        }
                    }
                }
            }
            Cut::Nothing | Cut::Unsupported | Cut::Rulings(_) => {}
        }
    }
    Some(out)
}

/// Whether record edge `i` lies in the plane (within its tolerance), and else where it crosses
/// it: its vertices on the plane, and the exact crossings between them.
fn edge_crossings(record: &Brep, i: usize, plane: &Plane, tiny: f64) -> (bool, Option<Vec<Crossing>>) {
    let Some(e) = record.edges.get(i) else { return (false, None) };
    let Ok(curve) = CurveEval::new(e.parameters()) else { return (false, None) };
    let [t0, t1] = e.interval;
    let vertex = |k: usize| e.vertices.get(k).and_then(|v| record.vertices.get(*v as usize)).map(|v| (Vec3::from(v.point), v.tolerance.max(tiny)));
    let (Some((p0, tol0)), Some((p1, tol1))) = (vertex(0), vertex(1)) else { return (false, None) };
    let etol = e.tolerance.max(tiny);
    let h = |t: f64| curve.point(t).map(|p| plane.height(p));
    let ts: Vec<f64> = (0..=SAMPLES).map(|k| t0 + (t1 - t0) * k as f64 / SAMPLES as f64).collect();
    let Some(hs) = ts.iter().map(|t| h(*t)).collect::<Option<Vec<f64>>>() else { return (false, None) };
    let (on0, on1) = (plane.height(p0).abs() <= tol0, plane.height(p1).abs() <= tol1);
    let inner = hs.get(1..SAMPLES).unwrap_or_default();
    if on0 && on1 && inner.iter().all(|x| x.abs() <= etol) {
        return (true, None);
    }
    let mut out = Vec::new();
    if on0 {
        out.push(Crossing { at: p0, tol: tol0 });
    }
    if on1 && p1.dist(p0) > tiny {
        out.push(Crossing { at: p1, tol: tol1 });
    }
    // Sign changes away from an end on the plane.
    for (k, (t, x)) in ts.windows(2).zip(hs.windows(2)).enumerate() {
        let (&[ta, tb], &[ha, hb]) = (t, x) else { continue };
        if (k == 0 && on0) || (k + 1 == SAMPLES && on1) {
            continue;
        }
        if ha == 0.0 && k > 0 {
            if let Some(p) = curve.point(ta) {
                out.push(Crossing { at: p, tol: etol });
            }
            continue;
        }
        if ha.signum() == hb.signum() || hb == 0.0 {
            continue;
        }
        let Some(t) = root(&h, ta, tb, ha) else { continue };
        if let Some(p) = curve.point(t) {
            out.push(Crossing { at: p, tol: etol });
        }
    }
    (false, Some(out))
}

/// The parameter in `[a, b]` where `h` is zero (it changes sign there), to the last bit.
fn root(h: &impl Fn(f64) -> Option<f64>, mut a: f64, mut b: f64, mut ha: f64) -> Option<f64> {
    for _ in 0..200 {
        let m = 0.5 * (a + b);
        if m <= a || m >= b {
            break;
        }
        let hm = h(m)?;
        if hm == 0.0 {
            return Some(m);
        }
        if hm.signum() == ha.signum() {
            (a, ha) = (m, hm);
        } else {
            b = m;
        }
    }
    // The end nearer the plane.
    let (x, y) = (h(a)?.abs(), h(b)?.abs());
    Some(if x <= y { a } else { b })
}

/// How far `p` is from edge `e` (sampled, then refined to the nearest point).
fn distance_to_edge(e: &solidcore_record::brep::Edge, p: Vec3) -> Option<f64> {
    let curve = CurveEval::new(e.parameters()).ok()?;
    let [t0, t1] = e.interval;
    let d = |t: f64| curve.point(t).map(|q| q.dist(p));
    let step = (t1 - t0) / SAMPLES as f64;
    let (mut best, mut bd) = (t0, f64::INFINITY);
    for k in 0..=SAMPLES {
        let t = t0 + step * k as f64;
        let x = d(t)?;
        if x < bd {
            (best, bd) = (t, x);
        }
    }
    // Golden section about the best sample.
    let (mut a, mut b) = ((best - step).max(t0), (best + step).min(t1));
    let g = 0.5 * (5f64.sqrt() - 1.0);
    for _ in 0..100 {
        let (c, e2) = (b - g * (b - a), a + g * (b - a));
        if d(c)? < d(e2)? {
            b = e2;
        } else {
            a = c;
        }
    }
    Some(bd.min(d(0.5 * (a + b))?))
}

/// What `plane` (unit normal `n`) cuts from surface `s`.
fn surface_cut(s: &Surface, plane: &Plane, n: Vec3, tiny: f64) -> Cut {
    let v = Vec3::from;
    let frame = |f: &solidcore_record::brep::Frame| {
        let (x, y) = (v(f.x_axis), v(f.y_axis));
        (v(f.origin), x.cross(y).normalized())
    };
    let local_dir = |d: Vec3| Vec2::new(d.dot(plane.x), d.dot(plane.y)).normalized();
    let flat = |origin: Vec3, m: Vec3| -> Cut {
        let Some(m) = m.normalized() else { return Cut::Unsupported };
        if m.cross(n).len() <= PARALLEL {
            return Cut::Nothing;
        }
        // Points o + a x + b y with m·(o + a x + b y − origin) = 0.
        let g = Vec2::new(m.dot(plane.x), m.dot(plane.y));
        let c = m.dot(origin - plane.origin);
        let Some(d) = g.perp().normalized() else { return Cut::Unsupported };
        Cut::Lines(vec![(g * (c / g.len2()), d)])
    };
    match s {
        Surface::Plane { origin, u_axis, v_axis } => flat(v(*origin), v(*u_axis).cross(v(*v_axis))),
        Surface::PlaneThrough { origin, u_point, v_point } => {
            let o = v(*origin);
            flat(o, (v(*u_point) - o).cross(v(*v_point) - o))
        }
        Surface::Cylinder { frame: f, radius } => {
            let (o, Some(z)) = frame(f) else { return Cut::Unsupported };
            if z.cross(n).len() <= PARALLEL {
                return Cut::Circles(vec![(plane.to_local(o), *radius)]);
            }
            if z.dot(n).abs() > PARALLEL {
                return Cut::Unsupported;
            }
            let h = plane.height(o);
            if h.abs() > *radius + tiny {
                return Cut::Nothing;
            }
            let w = (radius * radius - h * h).max(0.0).sqrt();
            let (Some(across), Some(d)) = (z.cross(n).normalized(), local_dir(z)) else { return Cut::Unsupported };
            let foot = o - n * h;
            let mut lines = vec![(plane.to_local(foot + across * w), d)];
            if w > tiny {
                lines.push((plane.to_local(foot - across * w), d));
            }
            Cut::Lines(lines)
        }
        Surface::Sphere { frame: f, radius } => {
            let (o, _) = frame(f);
            sphere(plane, o, *radius, tiny)
        }
        Surface::BallPatch { center, radius } => sphere(plane, v(*center), *radius, tiny),
        Surface::Cone { frame: f, radius_at_origin, half_angle } => {
            let (o, Some(z)) = frame(f) else { return Cut::Unsupported };
            if z.cross(n).len() > PARALLEL {
                return Cut::Unsupported;
            }
            let at = -plane.height(o) / z.dot(n);
            let r = radius_at_origin + at * half_angle.tan();
            if r <= tiny { Cut::Nothing } else { Cut::Circles(vec![(plane.to_local(o + z * at), r)]) }
        }
        Surface::Torus { frame: f, major_radius, minor_radius, .. } => {
            let (o, Some(z)) = frame(f) else { return Cut::Unsupported };
            if z.cross(n).len() > PARALLEL {
                return Cut::Unsupported;
            }
            let at = -plane.height(o) / z.dot(n);
            if at.abs() > *minor_radius + tiny {
                return Cut::Nothing;
            }
            let w = (minor_radius * minor_radius - at * at).max(0.0).sqrt();
            let c = plane.to_local(o + z * at);
            let mut circles = vec![(c, major_radius + w)];
            if w > tiny && major_radius - w > tiny {
                circles.push((c, major_radius - w));
            }
            Cut::Circles(circles)
        }
        Surface::Spline { u_knots, v_knots, points } => {
            let net: Vec<Vec<[f64; 4]>> = points.iter().map(|r| r.iter().map(|p| [p[0], p[1], p[2], 1.0]).collect()).collect();
            rulings(u_knots, v_knots, &net, n, local_dir, tiny)
        }
        Surface::RationalSpline { u_knots, v_knots, homogeneous } => rulings(u_knots, v_knots, homogeneous, n, local_dir, tiny),
        _ => Cut::Unsupported,
    }
}

/// The cut of a spline surface that is a curve extruded straight (degree 1 with two rows of
/// control points one fixed step apart, of equal weights) along a direction in the plane.
fn rulings(u_knots: &[f64], v_knots: &[f64], net: &[Vec<[f64; 4]>], n: Vec3, local_dir: impl Fn(Vec3) -> Option<Vec2>, tiny: f64) -> Cut {
    let point = |h: &[f64; 4]| Vec3::new(h[0] / h[3], h[1] / h[3], h[2] / h[3]);
    // The step between the two rows of each pair, when it is one step.
    let step = |pairs: &mut dyn Iterator<Item = (&[f64; 4], &[f64; 4])>| -> Option<Vec3> {
        let mut d: Option<Vec3> = None;
        for (a, b) in pairs {
            if (a[3] - b[3]).abs() > 1e-12 * a[3].abs() {
                return None;
            }
            let x = point(b) - point(a);
            if d.is_some_and(|d| d.dist(x) > tiny) {
                return None;
            }
            d = d.or(Some(x));
        }
        d.filter(|d| d.len() > tiny)
    };
    let along_v = (v_knots.len() == 4 && net.iter().all(|r| r.len() == 2)).then(|| step(&mut net.iter().map(|r| (&r[0], &r[1])))).flatten();
    let along_u = || (u_knots.len() == 4 && net.len() == 2).then(|| step(&mut net[0].iter().zip(&net[1]))).flatten();
    let Some(z) = along_v.or_else(along_u).and_then(|d| d.normalized()) else { return Cut::Unsupported };
    if z.dot(n).abs() > PARALLEL {
        return Cut::Unsupported;
    }
    local_dir(z).map_or(Cut::Unsupported, Cut::Rulings)
}

fn sphere(plane: &Plane, o: Vec3, radius: f64, tiny: f64) -> Cut {
    let h = plane.height(o);
    if h.abs() >= radius - tiny {
        return Cut::Nothing;
    }
    Cut::Circles(vec![(plane.to_local(o), (radius * radius - h * h).sqrt())])
}

/// Edge `e` of `body` seen along `plane`'s normal when it lies on a cylinder about that normal
/// and is no line or circle itself: the cylinder's arc between the projected ends of the edge,
/// on the side its middle lies (an edge with both ends at one point, the whole circle).
pub(super) fn cylinder_arc(body: &ModelBody, plane: &Plane, e: usize) -> Option<LinkGeom> {
    use solidcore_record::brep::Curve;
    if body.body.is_mesh() {
        return None;
    }
    let record = body.body.record();
    let edge = record.edges.get(e)?;
    if matches!(edge.curve, Curve::Line { .. } | Curve::Circle { .. }) {
        return None;
    }
    let n = plane.normal().normalized()?;
    let tiny = body.body.size().max(1e-9) * 1e-9;
    let (c, r) = edge.uses.iter().find_map(|u| match &record.faces.get(u.face as usize)?.surface {
        Surface::Cylinder { frame, radius } => {
            let z = Vec3::from(frame.x_axis).cross(Vec3::from(frame.y_axis)).normalized()?;
            (z.cross(n).len() <= PARALLEL).then(|| (plane.to_local(Vec3::from(frame.origin)), *radius))
        }
        _ => None,
    })?;
    let end = |k: usize| -> Option<(Vec2, f64)> {
        let v = record.vertices.get(*edge.vertices.get(k)? as usize)?;
        Some((plane.to_local(Vec3::from(v.point)), v.tolerance.max(tiny)))
    };
    let ((a, ta), (b, tb)) = (end(0)?, end(1)?);
    if (a.dist(c) - r).abs() > ta || (b.dist(c) - r).abs() > tb {
        return None;
    }
    if a.dist(b) <= ta.max(tb) {
        return Some(LinkGeom::Circle(c, r));
    }
    let [t0, t1] = edge.interval;
    let m = plane.to_local(CurveEval::new(edge.parameters()).ok()?.point(0.5 * (t0 + t1))?);
    let turn = |p: Vec2| ((p - c).angle() - (a - c).angle()).rem_euclid(std::f64::consts::TAU);
    Some(if turn(m) < turn(b) { LinkGeom::Arc { c, a, b } } else { LinkGeom::Arc { c, a: b, b: a } })
}

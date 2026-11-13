//! A circle or arc seen edge-on projects to a segment, exactly.

use super::*;

fn ring(c: Vec3, e1: Vec3, e2: Vec3, r: f64, from: f64, to: f64, n: usize) -> Vec<Vec3> {
    (0..=n).map(|i| from + (to - from) * i as f64 / n as f64).map(|t| c + (e1 * t.cos() + e2 * t.sin()) * r).collect()
}

fn line(g: &[LinkGeom]) -> (Vec2, Vec2) {
    match g {
        [LinkGeom::Line(a, b)] => (*a, *b),
        other => panic!("not one line: {other:?}"),
    }
}

fn same_ends((a, b): (Vec2, Vec2), p: Vec2, q: Vec2) -> bool {
    (a.dist(p) < 1e-12 && b.dist(q) < 1e-12) || (a.dist(q) < 1e-12 && b.dist(p) < 1e-12)
}

/// Rahmenrohr's PCB sketch: a tube's r 11 end circle in YZ, seen from the XY plane, is the
/// segment across its diameter to the last digit (a polyline fit put one end 1.4e-3 off).
#[test]
fn a_circle_seen_edge_on_is_its_diameter() {
    let plane = Plane { origin: Vec3::ZERO, x: Vec3::X, y: Vec3::Y };
    let c = Vec3::new(21.400572090709183, -106.06601717798213, 0.0);
    let pts = ring(c, Vec3::X, Vec3::Z, 11.0, 0.0, std::f64::consts::TAU, 96);
    let got = line(&fit_world(&plane, &pts));
    assert!(same_ends(got, Vec2::new(10.400572090709183, -106.06601717798213), Vec2::new(32.40057209070918, -106.06601717798213)), "{got:?}");
}

/// An arc seen edge-on: its ends and, where it passes them, the circle's extremes.
#[test]
fn an_arc_seen_edge_on_spans_its_extent() {
    let plane = Plane { origin: Vec3::ZERO, x: Vec3::X, y: Vec3::Y };
    // A half turn from (0,0,5) over (5,0,0) to (0,0,-5): its ends project onto the centre.
    let pts = ring(Vec3::ZERO, Vec3::Z, Vec3::X, 5.0, 0.0, std::f64::consts::PI, 64);
    assert!(same_ends(line(&fit_world(&plane, &pts)), Vec2::new(0.0, 0.0), Vec2::new(5.0, 0.0)));
    // 30 to 60 degrees off +Z towards +X: no extreme inside, its ends alone.
    let (a, b) = (30f64.to_radians(), 60f64.to_radians());
    let pts = ring(Vec3::ZERO, Vec3::Z, Vec3::X, 5.0, a, b, 32);
    assert!(same_ends(line(&fit_world(&plane, &pts)), Vec2::new(5.0 * a.sin(), 0.0), Vec2::new(5.0 * b.sin(), 0.0)));
}

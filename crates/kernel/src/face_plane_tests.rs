//! A face is a plane by its surface, not by its triangles: a curved side that bulges less than the
//! mesh tolerance meshes one triangle wide, all its triangles in one plane, and is still curved.

use solidcore_record::{Loop2, Plane, Region2, Seg2, Vec2};

/// A unit square whose bottom is a cubic and whose top is an arc, each bulging 3e-4 out of
/// line, extruded 10: four flat faces (the caps and the straight sides) and two curved ones.
#[test]
fn a_curve_flatter_than_the_mesh_tolerance_is_no_plane() {
    let (s, sag): (f64, f64) = (4e-4, 3e-4);
    let r = (0.25 + sag * sag) / (2.0 * sag);
    let half = (0.5 / r).asin();
    let segs = vec![
        Seg2::Cubic { p0: Vec2::new(0.0, 0.0), p1: Vec2::new(1.0 / 3.0, -s), p2: Vec2::new(2.0 / 3.0, -s), p3: Vec2::new(1.0, 0.0) },
        Seg2::Line { a: Vec2::new(1.0, 0.0), b: Vec2::new(1.0, 1.0) },
        Seg2::Arc { center: Vec2::new(0.5, 1.0 + sag - r), radius: r, start: std::f64::consts::FRAC_PI_2 - half, sweep: 2.0 * half },
        Seg2::Line { a: Vec2::new(0.0, 1.0), b: Vec2::new(0.0, 0.0) },
    ];
    let body = crate::extrude(&Plane::XY, &[Region2 { outer: Loop2 { segs }, holes: vec![] }], 0.0, 10.0).unwrap().remove(0);
    let faces = body.faces(1e-3).unwrap();
    let flat: Vec<_> = faces.iter().filter(|f| f.plane_normal.is_some()).map(|f| f.centroid).collect();
    assert_eq!(flat.len(), 4, "planes at {flat:?}");
    assert!(flat.iter().all(|c| (c.y - 0.5).abs() < 0.01), "a curved side taken for a plane: {flat:?}");
}

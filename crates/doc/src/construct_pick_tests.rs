//! How a reference's bare point (one without a name) picks its face.

use super::*;

/// The face the point lies on, not the first body's with a face within 5 % of its size: here a
/// rod's round side 0.5 mm from a pick on a box's side.
#[test]
fn a_bare_face_pick_is_the_face_the_point_lies_on() {
    let mut st = ModelState::default();
    st.bodies.push(ModelBody::new("Rod".into(), kernel::cylinder(Vec3::new(15.5, 5.0, 0.0), Vec3::Z, 5.0, 10.0).unwrap(), 0));
    st.bodies.push(ModelBody::new("Box".into(), kernel::box_solid(Vec3::ZERO, Vec3::new(10.0, 10.0, 10.0)).unwrap(), 0));
    let doc = Document::new("picks");
    let vals = BTreeMap::new();
    let ctx = Ctx { doc: &doc, vals: &vals, st: &st, depth: 0 };
    let p = plane(&ctx, &GeoRef::Face { at: Vec3::new(10.0, 5.0, 5.0), name: None }).unwrap();
    assert!(p.normal().x.abs() > 1.0 - 1e-9 && (p.origin - Vec3::new(10.0, 5.0, 5.0)).dot(p.normal()).abs() < 1e-9, "{p:?}");
    // A point on no face picks none, however near faces lie.
    let e = plane(&ctx, &GeoRef::Face { at: Vec3::new(10.3, 5.0, 5.0), name: None }).unwrap_err();
    assert!(e.to_string().contains("no face"), "{e}");
}

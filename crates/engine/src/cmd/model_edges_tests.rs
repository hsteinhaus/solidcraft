//! `model.edges` gives each edge's middle on the edge's own curve, not on a tessellation of it.

use serde_json::{Value, json};
use solidcore_record::Vec3;
use solidcraft_geom::prelude::*;

use crate::Session;

fn run(s: &mut Session, id: &str, p: Value) -> Value {
    match s.execute(id, &p) {
        Ok(v) => v,
        Err(e) => panic!("{id} {p}: {e}"),
    }
}

fn at(v: &Value) -> Vec3 {
    Vec3::new(v[0].as_f64().unwrap(), v[1].as_f64().unwrap(), v[2].as_f64().unwrap())
}

/// A 200 x 200 x 1 plate with one upright edge rounded at radius 1: the body is large enough
/// that the quarter arcs of the round (0.29 from their chords) tessellate as their chords.
#[test]
fn the_middle_of_a_small_arc_on_a_large_body_is_on_the_arc() {
    let mut s = Session::default();
    run(&mut s, "solid.box", json!({"length": 200, "width": 200, "height": 1}));
    let edges = run(&mut s, "model.edges", json!({"body": "Body1"}));
    let upright = edges["edges"].as_array().unwrap().iter().find(|e| (at(&e["start"]).z - at(&e["end"]).z).abs() > 0.5).unwrap();
    run(&mut s, "solid.fillet", json!({"edges": [upright["mid"]], "radius": 1}));
    let edges = run(&mut s, "model.edges", json!({"body": "Body1"}));
    let arcs: Vec<_> = edges["edges"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|e| (at(&e["start"]).dist(at(&e["end"])) - 2f64.sqrt()).abs() < 1e-9 && (at(&e["start"]).z - at(&e["end"]).z).abs() < 1e-9)
        .collect();
    assert_eq!(arcs.len(), 2, "the round's arcs on the top and the bottom: {edges}");
    for e in arcs {
        let chord = (at(&e["start"]) + at(&e["end"])) * 0.5;
        let rise = at(&e["mid"]).dist(chord);
        assert!((rise - (1.0 - 0.5f64.sqrt())).abs() < 1e-9, "the middle is {rise} from the chord: {e}");
    }
}

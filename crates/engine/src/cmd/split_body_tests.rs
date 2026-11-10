//! Split Body by a whole solid body (Fusion's Split Body with a body as its splitting tool, not
//! extended): each target splits into the pieces outside the tool and the pieces inside it; the
//! tool, here a body of a linked design's placement, stays as it is.

use serde_json::{Value, json};

use crate::Session;

fn run(s: &mut Session, id: &str, p: Value) -> Value {
    match s.execute(id, &p) {
        Ok(v) => v,
        Err(e) => panic!("{id} {p}: {e}"),
    }
}

/// Two 100 x 100 x 10 plates split by a plus of two 10 mm bars through both (a linked design's
/// body, placed 200 mm back along X): each plate into its four corners (45 x 45 x 10, 6 faces) and
/// the plus inside the tool (19000 mm³, 14 faces), corners first, the first keeping the plate's name.
#[test]
fn a_linked_placements_body_splits_plates_into_the_pieces_outside_and_inside_it() {
    let dir = std::env::temp_dir().join(format!("solidcraft-split-by-body-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("plus.solidcraft");
    let mut plus = Session::default();
    run(&mut plus, "solid.box", json!({"length": 100, "width": 10, "height": 40, "corner": [200, 45, -10], "body_name": "Plus"}));
    run(&mut plus, "solid.box", json!({"length": 10, "width": 100, "height": 40, "corner": [245, 0, -10], "operation": "join"}));
    run(&mut plus, "file.save", json!({"path": path}));
    let mut s = Session::default();
    run(&mut s, "solid.box", json!({"length": 100, "width": 100, "height": 10, "body_name": "Plate"}));
    run(&mut s, "solid.box", json!({"length": 100, "width": 100, "height": 10, "corner": [0, 0, 15], "body_name": "Slab", "operation": "new"}));
    let linked = run(&mut s, "component.insert_design", json!({"path": path, "translate": [-200, 0, 0]}));
    let _ = std::fs::remove_dir_all(&dir);
    let out = run(&mut s, "solid.split_body", json!({"bodies": ["Plate", "Slab"], "tool": {"body": "Plus", "occurrence": linked["occurrence"]}}));
    for target in ["Plate", "Slab"] {
        let pieces = out["pieces"][target].as_array().unwrap_or_else(|| panic!("{out}"));
        let got: Vec<(String, i64, u64)> = pieces
            .iter()
            .map(|p| (p["name"].as_str().unwrap().to_string(), p["volume"].as_f64().unwrap().round() as i64, p["faces"].as_u64().unwrap()))
            .collect();
        let name = |k: usize| if k == 0 { target.to_string() } else { format!("{target} ({k})") };
        let want: Vec<(String, i64, u64)> = (0..5).map(|k| (name(k), if k < 4 { 20250 } else { 19000 }, if k < 4 { 6 } else { 14 })).collect();
        assert_eq!(got, want, "{target}");
    }
    let tool = run(&mut s, "inspect.measure", json!({"bodies": ["Plus"]}));
    assert!((tool["bodies"][0]["volume_mm3"].as_f64().unwrap() - 76000.0).abs() < 1e-6, "the tool stays: {tool}");
}

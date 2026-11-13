//! Component patterns (`solid.pattern.*` with `components`): a Pattern feature that owns the
//! placements it makes ([`solidcraft_doc::PlacementPattern`]). Placements are the design's, not
//! the model's, so the session keeps them in step with the feature after each evaluation, before
//! the joints are solved: a place the pattern no longer has loses its placement (a joint on it
//! then fails), a new place gets one, and a placement still where the pattern last put it follows
//! the pattern. One a joint or a move put elsewhere stays there, as in Fusion, where the joints
//! after a pattern move its instances.

use solidcraft_doc::{FeatureKind, Mat, PatternInstance};

use crate::Session;

/// What [`sync`] does to one pattern source: the instances it keeps (updated), and the
/// occurrences to make (key, transform) and to delete.
struct Plan {
    feature: usize,
    source: usize,
    keep: Vec<PatternInstance>,
    make: Vec<([u32; 2], Mat)>,
    delete: Vec<u64>,
    moved: Vec<(u64, Mat)>,
}

/// Bring every live component pattern's placements in step with it. Whether anything changed.
pub(crate) fn sync(s: &mut Session) -> bool {
    let st = s.model.state();
    let mut plans = Vec::new();
    for (i, f) in s.doc.features.iter().enumerate() {
        let FeatureKind::Pattern { pattern, placements, .. } = &f.kind else { continue };
        // A rolled back, suppressed or failing pattern leaves its placements as they are.
        let live = s.doc.marker.is_none_or(|m| i < m) && !f.suppressed && s.model.result(f.id).is_some_and(|r| r.error.is_none() && !r.skipped);
        if placements.is_empty() || !live {
            continue;
        }
        let Ok(mats) = solidcraft_doc::pattern_matrices(&s.doc, &st, pattern) else { continue };
        let Ok(keys) = solidcraft_doc::pattern_keys(&s.doc, pattern, mats.len()) else { continue };
        let wanted: std::collections::BTreeMap<[u32; 2], Mat> = keys.into_iter().zip(mats).collect();
        for (k, pp) in placements.iter().enumerate() {
            let Some(src) = s.doc.occurrences.iter().find(|o| o.id == pp.source) else { continue };
            let at = |m: &Mat| solidcraft_doc::mat_mul(m, &src.transform);
            let mut plan = Plan { feature: i, source: k, keep: Vec::new(), make: Vec::new(), delete: Vec::new(), moved: Vec::new() };
            for inst in &pp.instances {
                let occ = s.doc.occurrences.iter().find(|o| o.id == inst.occurrence);
                match (wanted.get(&inst.key), occ) {
                    (None, Some(_)) => plan.delete.push(inst.occurrence),
                    (None, None) => {}
                    (Some(m), Some(o)) if o.transform == inst.placed && at(m) != inst.placed => {
                        plan.moved.push((o.id, at(m)));
                        plan.keep.push(PatternInstance { placed: at(m), ..inst.clone() });
                    }
                    // Deleted (not made again), or kept where it is.
                    (Some(_), _) => plan.keep.push(inst.clone()),
                }
            }
            plan.make = wanted.iter().filter(|(key, _)| !pp.instances.iter().any(|x| x.key == **key)).map(|(key, m)| (*key, at(m))).collect();
            if !plan.make.is_empty() || !plan.delete.is_empty() || !plan.moved.is_empty() {
                plans.push(plan);
            }
        }
    }
    if plans.is_empty() {
        return false;
    }
    let doc = s.doc_mut();
    for plan in plans {
        doc.occurrences.retain(|o| !plan.delete.contains(&o.id));
        for (id, m) in &plan.moved {
            if let Some(o) = doc.occurrences.iter_mut().find(|o| o.id == *id) {
                o.transform = *m;
            }
        }
        let mut keep = plan.keep;
        let Some(FeatureKind::Pattern { placements, .. }) = doc.features.get(plan.feature).map(|f| &f.kind) else { continue };
        let Some(source) = placements.get(plan.source).and_then(|pp| doc.occurrences.iter().find(|o| o.id == pp.source)).cloned() else { continue };
        for (key, m) in plan.make {
            // Refused (too many occurrences): the place stays empty until the next sync.
            if let Ok(occurrence) = doc.add_occurrence(source.component, source.parent, m) {
                keep.push(PatternInstance { key, occurrence, placed: m });
            }
        }
        keep.sort_by_key(|x| x.key);
        if let Some(FeatureKind::Pattern { placements, .. }) = doc.features.get_mut(plan.feature).map(|f| &mut f.kind)
            && let Some(pp) = placements.get_mut(plan.source)
        {
            pp.instances = keep;
        }
    }
    true
}

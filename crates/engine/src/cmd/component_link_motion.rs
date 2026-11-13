//! Motion of a linked design: the parts its joints move, each in a component of its own placed
//! beside every placement of the linked component, joined to it by this design's own copies of
//! the linked design's joints (on the linked faces, by name), so each placement moves on its own.

use std::collections::{BTreeMap, BTreeSet};

use solidcraft_doc::joints::{JointKind, JointOrigin, Snap};
use solidcraft_doc::{DesignLink, LinkedEnd, LinkedJoint, LinkedPart, LinkedPlacement};

use crate::{Result, Session};

/// What moves in `part`, a design being linked: its parts (by the linked bodies' indices) and
/// its joints between them, and the part each occurrence moves with (0: the base). `sources`
/// are the world names of the bodies the linked bodies come from, by index.
pub(super) fn design_motion(
    part: &Session,
    sources: &[String],
    warnings: &mut Vec<String>,
) -> (Vec<LinkedPart>, Vec<LinkedJoint>, BTreeMap<u64, usize>) {
    let doc = &part.doc;
    let st = part.model.state();
    let joints: Vec<_> = doc.assembly.joints.iter().filter(|j| !j.suppressed).collect();
    let moving = |k: JointKind| k != JointKind::Rigid;
    if !joints.iter().any(|j| moving(j.kind)) {
        return Default::default();
    }
    // Rigid sets: occurrences joined rigidly, and an occurrence no joint places with the one it
    // sits in.
    let mut up: BTreeMap<u64, u64> = BTreeMap::new();
    fn find(up: &mut BTreeMap<u64, u64>, x: u64) -> u64 {
        let mut r = x;
        while let Some(&p) = up.get(&r).filter(|p| **p != r) {
            r = p;
        }
        up.insert(x, r);
        r
    }
    let mut join = |up: &mut BTreeMap<u64, u64>, a: u64, b: u64| {
        let (ra, rb) = (find(up, a), find(up, b));
        if ra != rb {
            up.insert(rb, ra);
        }
    };
    let jointed: BTreeSet<u64> = joints.iter().flat_map(|j| [j.a.occurrence, j.b.occurrence]).collect();
    for o in &doc.occurrences {
        find(&mut up, o.id);
        if !jointed.contains(&o.id) {
            let parent = if o.parent == 0 { 0 } else { doc.occurrence_of(o.parent).map_or(0, |p| p.id) };
            join(&mut up, parent, o.id);
        }
    }
    for j in joints.iter().filter(|j| !moving(j.kind)) {
        join(&mut up, j.a.occurrence, j.b.occurrence);
    }
    // The base: what stays put (the root, grounded occurrences, else the solve's first fixed).
    let base = find(&mut up, 0);
    for f in solidcraft_doc::joints::solve(doc, &st).fixed {
        join(&mut up, base, f);
    }
    let base = find(&mut up, 0);
    // Parts, in the order the moving joints reach them from the base.
    let mut part_of: BTreeMap<u64, usize> = BTreeMap::from([(base, 0)]);
    let mut order = vec![base];
    let mut grew = true;
    while grew {
        grew = false;
        for j in joints.iter().filter(|j| moving(j.kind)) {
            let (a, b) = (find(&mut up, j.a.occurrence), find(&mut up, j.b.occurrence));
            for (from, to) in [(a, b), (b, a)] {
                if part_of.contains_key(&from) && !part_of.contains_key(&to) {
                    part_of.insert(to, order.len());
                    order.push(to);
                    grew = true;
                }
            }
        }
    }
    if order.len() < 2 {
        return Default::default();
    }
    let occurrence_part = |up: &mut BTreeMap<u64, u64>, occ: u64| {
        let r = find(up, occ);
        part_of.get(&r).copied().unwrap_or(0)
    };
    // Each part named after its shallowest occurrence's component.
    let depth = |occ: u64| doc.first_path(doc.occurrences.iter().find(|o| o.id == occ).map_or(0, |o| o.parent)).len();
    let mut parts: Vec<LinkedPart> = (1..order.len()).map(|_| LinkedPart::default()).collect();
    let mut shallowest: Vec<Option<usize>> = vec![None; parts.len()];
    let mut by_occ = BTreeMap::new();
    for o in &doc.occurrences {
        let k = occurrence_part(&mut up, o.id);
        by_occ.insert(o.id, k);
        if let Some(slot) = k.checked_sub(1).and_then(|i| parts.get_mut(i).zip(shallowest.get_mut(i))) {
            let d = depth(o.id);
            if slot.1.is_none_or(|was| d < was) {
                *slot.1 = Some(d);
                slot.0.name = doc.components.iter().find(|c| c.id == o.component).map_or_else(|| o.name.clone(), |c| c.name.clone());
            }
        }
    }
    // Told apart by name, which keeps each part's component when the link is read again.
    for i in 1..parts.len() {
        let taken = |n: &str, parts: &[LinkedPart]| parts.iter().take(i).any(|p| p.name == n);
        let stem = parts.get(i).map(|p| p.name.clone()).unwrap_or_default();
        let mut name = stem.clone();
        for k in 2.. {
            if !taken(&name, &parts) {
                break;
            }
            name = format!("{stem} ({k})");
        }
        if let Some(p) = parts.get_mut(i) {
            p.name = name;
        }
    }
    // Bodies by where they are placed.
    let mut body_part: BTreeMap<String, usize> = BTreeMap::new();
    for (name, _, path, _) in doc.body_frames(&st) {
        body_part.insert(name, path.last().map_or(0, |o| by_occ.get(o).copied().unwrap_or(0)));
    }
    for (k, source) in sources.iter().enumerate() {
        if let Some(p) = body_part.get(source).and_then(|p| p.checked_sub(1)).and_then(|i| parts.get_mut(i)) {
            p.bodies.push(k);
        }
    }
    // A linked body by the world name of the body a joint origin is on.
    let linked_body = |o: &JointOrigin| -> Option<usize> {
        let body = o.body.as_deref()?;
        let comp = doc.occurrences.iter().find(|x| x.id == o.occurrence).map_or(0, |x| x.component);
        let mut path = doc.first_path(doc.occurrences.iter().find(|x| x.id == o.occurrence).map_or(0, |x| x.parent));
        if o.occurrence != 0 {
            path.push(o.occurrence);
        }
        let world = doc.world_body_name(&st, body, comp, &path);
        sources.iter().position(|s| *s == world)
    };
    let end = |up: &mut BTreeMap<u64, u64>, o: &JointOrigin| -> Option<LinkedEnd> {
        let m = if o.occurrence == 0 { solidcraft_doc::IDENTITY } else { solidcraft_doc::joints::occurrence_world(doc, o.occurrence) };
        Some(LinkedEnd { part: occurrence_part(up, o.occurrence), body: linked_body(o), snap: o.snap.mapped(&m)? })
    };
    let mut linked = Vec::new();
    for j in joints.iter().filter(|j| moving(j.kind)) {
        let (Some(a), Some(b)) = (end(&mut up, &j.a), end(&mut up, &j.b)) else {
            warnings.push(format!("joint {} does not move the linked design: one of its origins is not on its bodies' faces", j.name));
            continue;
        };
        if a.part == b.part {
            continue;
        }
        linked.push(LinkedJoint {
            name: j.name.clone(),
            kind: j.kind,
            a,
            b,
            values: j.values.clone(),
            limits: j.limits.clone(),
            flip: j.flip,
            axis: j.axis,
            offset: j.offset,
            angle: j.angle,
        });
    }
    (parts, linked, by_occ)
}

/// Bring this design in line with the moving parts the link of `comp` keeps: the parts'
/// components hold their bodies, each placement of `comp` has the parts placed beside it and the
/// linked joints between them, and joints and origins on a placement go to the part they are on.
pub(super) fn sync(s: &mut Session, comp: u64) -> Result<()> {
    let Some(mut link) = s.doc.links.get(&comp).cloned() else { return Ok(()) };
    if link.parts.is_empty() {
        if !link.placements.is_empty() {
            let gone: BTreeSet<u64> = link.placements.iter().flat_map(|p| p.parts.iter().copied()).collect();
            let joints: BTreeSet<u64> = link.placements.iter().flat_map(|p| p.joints.iter().copied()).collect();
            let doc = s.doc_mut();
            doc.assembly.joints.retain(|j| !joints.contains(&j.id));
            doc.remove_occurrences(&gone);
            if let Some(l) = doc.links.get_mut(&comp) {
                l.placements.clear();
            }
            s.refresh();
        }
        return Ok(());
    }
    let st = s.model.state();
    let names = body_names(s, &link, &st);
    let doc = s.doc_mut();
    let cname = doc.components.iter().find(|c| c.id == comp).map(|c| (c.name.clone(), c.parent));
    let Some((cname, cparent)) = cname else { return Ok(()) };
    // The parts' components, and their bodies.
    let old: BTreeSet<u64> = doc.body_components.iter().filter(|(_, c)| link.parts.iter().any(|p| p.component == **c)).map(|(_, c)| *c).collect();
    doc.body_components.retain(|_, c| !old.contains(c));
    for part in &mut link.parts {
        if part.component == 0 || !doc.components.iter().any(|c| c.id == part.component) {
            part.component = doc.add_component(Some(&format!("{cname} {}", part.name)), cparent)?;
            // Placed only beside the link's placements, not where the new component is put.
            let own: BTreeSet<u64> = doc.occurrences.iter().filter(|o| o.component == part.component).map(|o| o.id).collect();
            doc.remove_occurrences(&own);
        }
        for name in part.bodies.iter().filter_map(|k| names.get(*k).cloned().flatten()) {
            doc.body_components.insert(name, part.component);
        }
    }
    // Placements.
    let placed: Vec<(u64, u64, solidcraft_doc::Mat, String)> =
        doc.occurrences.iter().filter(|o| o.component == comp).map(|o| (o.id, o.parent, o.transform, o.name.clone())).collect();
    let mut gone: BTreeSet<u64> = BTreeSet::new();
    let mut dropped_joints: BTreeSet<u64> = BTreeSet::new();
    let mut placements = Vec::new();
    for lp in std::mem::take(&mut link.placements) {
        if placed.iter().any(|(id, ..)| *id == lp.occurrence) {
            placements.push(lp);
        } else {
            gone.extend(lp.parts.iter().copied());
            dropped_joints.extend(lp.joints.iter().copied());
        }
    }
    for (occ, parent, transform, oname) in &placed {
        let at = match placements.iter().position(|p| p.occurrence == *occ) {
            Some(i) => i,
            None => {
                placements.push(LinkedPlacement { occurrence: *occ, ..Default::default() });
                placements.len() - 1
            }
        };
        let Some(lp) = placements.get_mut(at) else { continue };
        // A part's occurrence: kept while it still places that part, else made where the
        // placement is (the linked bodies' rest pose).
        let mut occs = Vec::new();
        for (i, part) in link.parts.iter().enumerate() {
            let kept = lp.parts.get(i).copied().filter(|q| doc.occurrences.iter().any(|o| o.id == *q && o.component == part.component));
            occs.push(match kept {
                Some(q) => q,
                None => doc.add_occurrence(part.component, *parent, *transform)?,
            });
        }
        gone.extend(lp.parts.iter().filter(|q| !occs.contains(q)));
        lp.parts = occs;
        let side = |e: &LinkedEnd, lp: &LinkedPlacement| JointOrigin {
            occurrence: if e.part == 0 { *occ } else { lp.parts.get(e.part - 1).copied().unwrap_or(*occ) },
            body: e.body.and_then(|k| names.get(k).cloned().flatten()),
            snap: linked_snap(&e.snap, link.feature),
        };
        let mut ids = Vec::new();
        for (i, lj) in link.joints.iter().enumerate() {
            let (a, b) = (side(&lj.a, lp), side(&lj.b, lp));
            let name = format!("{oname}/{}", lj.name);
            let existing = lp.joints.get(i).copied().and_then(|id| doc.assembly.joints.iter_mut().find(|j| j.id == id && j.kind == lj.kind));
            let id = match existing {
                Some(j) => {
                    (j.name, j.a, j.b, j.limits, j.flip, j.axis, j.offset, j.angle) =
                        (name, a, b, lj.limits.clone(), lj.flip, lj.axis, lj.offset, lj.angle);
                    if j.values.len() != lj.values.len() {
                        j.values.clone_from(&lj.values);
                    }
                    j.id
                }
                None => {
                    let id = doc.allocate_id()?;
                    doc.assembly.joints.push(solidcraft_doc::joints::Joint {
                        id,
                        name,
                        kind: lj.kind,
                        a,
                        b,
                        values: lj.values.clone(),
                        limits: lj.limits.clone(),
                        flip: lj.flip,
                        axis: lj.axis,
                        offset: lj.offset,
                        angle: lj.angle,
                        suppressed: false,
                        exprs: Default::default(),
                    });
                    id
                }
            };
            ids.push(id);
        }
        dropped_joints.extend(lp.joints.iter().filter(|id| !ids.contains(id)));
        lp.joints = ids;
    }
    doc.assembly.joints.retain(|j| !dropped_joints.contains(&j.id));
    doc.remove_occurrences(&gone);
    link.placements = placements;
    // Joints and origins on a placement go with the part their origin or body is on.
    let link_ref = &link;
    let generated: BTreeSet<u64> = link.placements.iter().flat_map(|p| p.joints.iter().copied()).collect();
    let mut asm = doc.assembly.clone();
    let mut moved = false;
    for j in asm.joints.iter_mut().filter(|j| !generated.contains(&j.id)) {
        for o in [&mut j.a, &mut j.b] {
            moved |= retarget(link_ref, &names, o);
        }
    }
    for o in &mut asm.origins {
        moved |= retarget(link_ref, &names, &mut o.origin);
    }
    if moved {
        doc.assembly = asm;
    }
    doc.links.insert(comp, link);
    s.refresh();
    let _ = super::super::joints::place(s);
    Ok(())
}

/// The link of `comp` takes the parts and joints its design has now (`parts`, `joints`), each
/// part keeping its component where the design still has a part of that name. A part it no longer
/// has goes with its component; what was joined to it goes back to its placement. Whether they
/// changed.
pub(super) fn follow(s: &mut Session, comp: u64, mut parts: Vec<LinkedPart>, joints: Vec<LinkedJoint>) -> Result<bool> {
    let Some(link) = s.doc.links.get(&comp) else { return Ok(false) };
    for p in &mut parts {
        p.component = link.parts.iter().find(|o| o.name == p.name).map_or(0, |o| o.component);
    }
    if link.parts == parts && link.joints == joints {
        return Ok(false);
    }
    let stale: Vec<u64> = link.parts.iter().map(|p| p.component).filter(|c| *c != 0 && !parts.iter().any(|p| p.component == *c)).collect();
    // Origins on a stale part's occurrence go back to the placement beside it.
    let back: BTreeMap<u64, u64> = link
        .placements
        .iter()
        .flat_map(|lp| lp.parts.iter().map(move |q| (*q, lp.occurrence)))
        .filter(|(q, _)| s.doc.occurrences.iter().any(|o| o.id == *q && stale.contains(&o.component)))
        .collect();
    let doc = s.doc_mut();
    if !back.is_empty() {
        let mut asm = doc.assembly.clone();
        let ends = asm.joints.iter_mut().flat_map(|j| [&mut j.a, &mut j.b]).chain(asm.origins.iter_mut().map(|o| &mut o.origin));
        for o in ends {
            if let Some(p) = back.get(&o.occurrence) {
                o.occurrence = *p;
            }
        }
        doc.assembly = asm;
    }
    for c in stale {
        doc.remove_component_tree(c)?;
    }
    if let Some(link) = doc.links.get_mut(&comp) {
        link.parts = parts;
        link.joints = joints;
    }
    Ok(true)
}

/// Move `o`, a joint origin on a placement of `link`'s component or of one of its parts, to the
/// occurrence of the part its linked origin or body moves with there. Whether it moved.
fn retarget(link: &DesignLink, names: &[Option<String>], o: &mut JointOrigin) -> bool {
    let Some(lp) = link.placements.iter().find(|p| p.occurrence == o.occurrence || p.parts.contains(&o.occurrence)) else { return false };
    let part = match (&o.snap, o.body.as_deref()) {
        (Snap::Linked { origin }, _) => link.origins.iter().find(|l| l.name == *origin).map(|l| l.part),
        (_, Some(body)) => {
            let k = names.iter().position(|n| n.as_deref() == Some(body));
            Some(k.and_then(|k| link.parts.iter().position(|p| p.bodies.contains(&k))).map_or(0, |i| i + 1))
        }
        _ => None,
    };
    let Some(part) = part else { return false };
    let want = if part == 0 { lp.occurrence } else { lp.parts.get(part - 1).copied().unwrap_or(lp.occurrence) };
    if want == o.occurrence {
        return false;
    }
    o.occurrence = want;
    true
}

/// The occurrence a linked origin `origin` given on `occ` is on: the part's beside it when the
/// origin moves with one.
pub(crate) fn origin_occurrence(doc: &solidcraft_doc::Document, occ: u64, origin: &str) -> u64 {
    let comp = doc.occurrences.iter().find(|o| o.id == occ).map_or(0, |o| o.component);
    let Some(link) = doc.links.get(&comp) else { return occ };
    let mut o = JointOrigin { occurrence: occ, body: None, snap: Snap::Linked { origin: origin.to_string() } };
    retarget(link, &[], &mut o);
    o.occurrence
}

/// `snap`, on the linked design's faces, as this design names them (the Import feature `fid`).
fn linked_snap(snap: &Snap, fid: u64) -> Snap {
    let mut out = snap.clone();
    for name in out.face_names_mut() {
        *name = solidcraft_doc::naming::link_face_name(fid, name);
    }
    out
}

/// This design's names of the linked bodies, by their index in the linked design's bodies
/// (`None` for one the Import feature doesn't build).
fn body_names(s: &Session, link: &DesignLink, st: &solidcraft_doc::ModelState) -> Vec<Option<String>> {
    let Some(solidcraft_doc::FeatureKind::Import { skip, .. }) = s.doc.feature(link.feature).map(|f| &f.kind) else { return Vec::new() };
    let built: Vec<String> = st.bodies.iter().filter(|b| b.feature == link.feature).map(|b| b.name.clone()).collect();
    let mut names = built.into_iter();
    let count = names.len() + skip.len();
    (0..count).map(|k| if skip.contains(&(k as u64)) { None } else { names.next() }).collect()
}

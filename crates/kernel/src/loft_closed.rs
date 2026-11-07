//! A surface lofted all the way round through open sections: the sections close up into a ring
//! (the last one runs on into the first), so the surface is periodic across them and smooth
//! everywhere, also where it is cut into faces.
//!
//! Each section must be one smooth curve (a spline's pieces are joined first, as the open loft
//! joins them). The sections are made compatible (one degree, one knot vector, by exact knot
//! insertion), and every control point is carried round the ring by periodic cubic interpolation
//! through the sections at equally spaced parameters: the surface passes through every section
//! exactly and is C² across all of them. Sections that are mirror images of each other about one
//! of them (a dome's radial profiles) give a surface that is mirror symmetric about that one, so
//! it leaves it square.
//!
//! The ring is cut into faces at the `seams` (section indices, in loft order): each face is the
//! surface between two seams, exactly (a chain of Bézier spans across the sections), and the faces
//! meet along the seam sections themselves. A boundary that collapses to a point (sections that all
//! end at one apex) becomes that point.

use std::collections::BTreeSet;

use solidcore_record::{Plane, Seg2};
use solidcraft_geom::bspline::{CurveN, HPoint, NurbsSurface};
use solidcraft_geom::prelude::*;
use solidcraft_geom::spline_space::{Space, embed};

use crate::body::Body;
use crate::{KernelError, Result, record};

fn invalid(what: impl Into<String>) -> KernelError {
    KernelError::Invalid(what.into())
}

fn failed(what: impl Into<String>) -> KernelError {
    KernelError::Failed(what.into())
}

/// Most sections a closed loft takes.
const MAX_SECTIONS: usize = 50;

/// A surface lofted through `sections` (each one open chain of sketch segments on its plane, all
/// running the same way) and back into the first, cut into faces at the section indices `seams`.
/// Walls only; the history has no operands: every face is new.
pub fn loft_sheet_closed(sections: &[(Plane, Vec<Seg2>)], seams: &[usize]) -> Result<Body> {
    let n = sections.len();
    if !(3..=MAX_SECTIONS).contains(&n) {
        return Err(invalid(format!("a closed loft needs 3…{MAX_SECTIONS} sections")));
    }
    let seams: BTreeSet<usize> = seams.iter().copied().collect();
    if seams.is_empty() || seams.iter().any(|s| *s >= n) {
        return Err(invalid("a closed loft's seams are section numbers, at least one, each less than the section count"));
    }
    for (plane, chain) in sections {
        if !crate::build::frame_ok(plane) {
            return Err(invalid("a loft section's frame is not finite and orthonormal"));
        }
        crate::surfaces::check_chains(std::slice::from_ref(chain), "loft")?;
        let closed = match (chain.first(), chain.last()) {
            (Some(f), Some(l)) => l.end().dist(f.start()) <= 1e-6,
            _ => return Err(invalid("a loft section has no curves")),
        };
        if closed {
            return Err(invalid("a closed loft runs through open sections (a closed profile lofts round as a tube)"));
        }
    }
    let matched = crate::surfaces::smooth_runs_joined(crate::surfaces::open_sections(sections)?);
    let curves: Vec<&CurveN<4>> = matched
        .iter()
        .map(|s| match s.curves.as_slice() {
            [c] => Ok(c),
            _ => Err(invalid("not supported yet: a closed loft through sections that are not each one smooth curve")),
        })
        .collect::<Result<_>>()?;
    let columns = compatible(&curves)?;
    let (u_degree, u_knots, rows) = columns;
    let m = rows.first().map_or(0, Vec::len);
    let ring = periodic(&rows)?;
    let seams: Vec<usize> = seams.into_iter().collect();
    let mut pieces = Vec::with_capacity(seams.len());
    for (k, &a) in seams.iter().enumerate() {
        let b = seams.get(k + 1).copied().unwrap_or_else(|| seams.first().copied().unwrap_or(a) + n);
        let span = if b > a { b - a } else { n };
        let net = patch(&ring, &rows, a, span, u_degree, &u_knots, m)?;
        let sheet = crate::surfaces::tensor_sheet(&net)?;
        pieces.push((sheet.record().clone(), vec![record::generated(Vec::new()); 1]));
    }
    crate::surfaces::sheet_of(pieces, true, &[])
}

/// The sections' curves on one degree and knot vector, each parameterised over `[0, 1]`: the
/// degree, the knots and each section's control points (homogeneous, unit weights).
fn compatible(curves: &[&CurveN<4>]) -> Result<(usize, Vec<f64>, Vec<Vec<HPoint>>)> {
    let first = curves.first().ok_or_else(|| invalid("no sections"))?;
    let p = first.degree;
    let polynomial = |c: &CurveN<4>| c.ctrl.iter().all(|h| (h[3] - 1.0).abs() <= 1e-12);
    if curves.iter().any(|c| c.degree != p || !polynomial(c)) {
        return Err(invalid("not supported yet: a closed loft through sections of different degrees or with weights"));
    }
    // Each curve's knots mapped onto [0, 1].
    let unit: Vec<Vec<f64>> = curves
        .iter()
        .map(|c| {
            let (t0, t1) = (*c.knots.get(p)?, *c.knots.get(c.ctrl.len())?);
            (t1 > t0).then(|| c.knots.iter().map(|t| ((t - t0) / (t1 - t0)).clamp(0.0, 1.0)).collect())
        })
        .collect::<Option<_>>()
        .ok_or_else(|| invalid("a section curve has an empty parameter range"))?;
    // Union of the interior knots, each at the highest multiplicity any curve has it.
    let mut merged: Vec<(f64, usize)> = Vec::new();
    for k in &unit {
        let interior = k.get(p + 1..k.len().saturating_sub(p + 1)).unwrap_or(&[]);
        let mut i = 0;
        while let Some(&t) = interior.get(i) {
            let mult = interior.iter().skip(i).take_while(|x| (**x - t).abs() <= 1e-10).count();
            match merged.iter_mut().find(|(x, _)| (*x - t).abs() <= 1e-10) {
                Some(e) => e.1 = e.1.max(mult),
                None => merged.push((t, mult)),
            }
            i += mult.max(1);
        }
    }
    merged.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut knots = vec![0.0; p + 1];
    for (t, mult) in &merged {
        knots.extend(std::iter::repeat_n(*t, *mult));
    }
    knots.extend(std::iter::repeat_n(1.0, p + 1));
    let target = Space::new(p, knots.clone()).map_err(|e| invalid(format!("closed loft: {e}")))?;
    let rows = curves
        .iter()
        .zip(&unit)
        .map(|(c, k)| {
            let from = Space::new(p, k.clone()).map_err(|e| invalid(format!("closed loft: {e}")))?;
            let ctrl: Vec<HPoint> = c.ctrl.iter().map(|h| HPoint { x: h[0], y: h[1], z: h[2], w: h[3] }).collect();
            embed(&from, &ctrl, &target).map_err(|e| failed(format!("closed loft: making the sections compatible: {e}")))
        })
        .collect::<Result<Vec<_>>>()?;
    Ok((p, knots, rows))
}

/// The periodic cubic B-spline control points (uniform knots, one per section) that pass through
/// every section's control points: for each section k, (X[k−1] + 4 X[k] + X[k+1]) / 6 = Q[k],
/// indices round the ring. Solved once for the circulant matrix, applied to every column.
fn periodic(rows: &[Vec<HPoint>]) -> Result<Vec<Vec<HPoint>>> {
    let n = rows.len();
    // The inverse of the circulant (1, 4, 1) / 6: its first row, by Gaussian elimination on the
    // n × n system (n is at most MAX_SECTIONS).
    let mut a = vec![0.0; n * n];
    for k in 0..n {
        for (d, w) in [(n - 1, 1.0 / 6.0), (0, 4.0 / 6.0), (1, 1.0 / 6.0)] {
            if let Some(x) = a.get_mut(k * n + (k + d) % n) {
                *x += w;
            }
        }
    }
    let inv = invert(n, a).ok_or_else(|| failed("closed loft: the ring of sections has no interpolant"))?;
    let m = rows.first().map_or(0, Vec::len);
    let mut out = vec![vec![HPoint::ZERO; m]; n];
    for (k, out_k) in out.iter_mut().enumerate() {
        for (i, x) in out_k.iter_mut().enumerate() {
            let mut s = HPoint::ZERO;
            for (j, row) in rows.iter().enumerate() {
                let w = inv.get(k * n + j).copied().unwrap_or(0.0);
                s += *row.get(i).ok_or_else(|| failed("closed loft: section sizes differ"))? * w;
            }
            *x = s;
        }
    }
    Ok(out)
}

/// The inverse of the n × n matrix `a` (row major), by Gauss-Jordan elimination with partial
/// pivoting; `None` when it is singular.
fn invert(n: usize, mut a: Vec<f64>) -> Option<Vec<f64>> {
    let mut inv = vec![0.0; n * n];
    for k in 0..n {
        *inv.get_mut(k * n + k)? = 1.0;
    }
    for c in 0..n {
        let size = |r: usize| a.get(r * n + c).map_or(0.0, |v| v.abs());
        let pivot = (c..n).max_by(|x, y| size(*x).total_cmp(&size(*y)))?;
        if a.get(pivot * n + c)?.abs() < 1e-12 {
            return None;
        }
        for j in 0..n {
            a.swap(c * n + j, pivot * n + j);
            inv.swap(c * n + j, pivot * n + j);
        }
        let d = *a.get(c * n + c)?;
        for j in 0..n {
            *a.get_mut(c * n + j)? /= d;
            *inv.get_mut(c * n + j)? /= d;
        }
        for r in 0..n {
            if r == c {
                continue;
            }
            let f = *a.get(r * n + c)?;
            if f == 0.0 {
                continue;
            }
            for j in 0..n {
                let (ac, ic) = (*a.get(c * n + j)?, *inv.get(c * n + j)?);
                *a.get_mut(r * n + j)? -= f * ac;
                *inv.get_mut(r * n + j)? -= f * ic;
            }
        }
    }
    Some(inv)
}

/// The surface between section `a` and the section `span` further round: u along the sections'
/// curves, v round the ring as a chain of Bézier spans (each span of the uniform periodic spline,
/// exactly), its first and last columns the two sections' own control points.
fn patch(ring: &[Vec<HPoint>], rows: &[Vec<HPoint>], a: usize, span: usize, p: usize, u_knots: &[f64], m: usize) -> Result<NurbsSurface> {
    let n = ring.len();
    let at = |k: usize, i: usize| ring.get(k % n).and_then(|r| r.get(i)).copied().ok_or_else(|| failed("closed loft: ring index"));
    let section = |k: usize, i: usize| rows.get(k % n).and_then(|r| r.get(i)).copied().ok_or_else(|| failed("closed loft: section index"));
    // v control columns: the first section, then each span's two inner Bézier points and its end
    // (the next section).
    let mut columns: Vec<Vec<HPoint>> = Vec::with_capacity(3 * span + 1);
    columns.push((0..m).map(|i| section(a, i)).collect::<Result<_>>()?);
    for s in 0..span {
        let k = a + s + n;
        let mut b1 = Vec::with_capacity(m);
        let mut b2 = Vec::with_capacity(m);
        for i in 0..m {
            let (x0, x1) = (at(k, i)?, at(k + 1, i)?);
            b1.push(x0 * (4.0 / 6.0) + x1 * (2.0 / 6.0));
            b2.push(x0 * (2.0 / 6.0) + x1 * (4.0 / 6.0));
        }
        columns.push(b1);
        columns.push(b2);
        columns.push((0..m).map(|i| section(k + 1, i)).collect::<Result<_>>()?);
    }
    let mut v_knots = vec![0.0; 4];
    for s in 1..span {
        v_knots.extend(std::iter::repeat_n(s as f64, 3));
    }
    v_knots.extend(std::iter::repeat_n(span as f64, 4));
    // The sections are polynomial (`compatible`), so is the surface: its weights are exactly 1,
    // not the sums of the interpolation's coefficients, which round off it.
    let mut ctrl = Vec::with_capacity(m * columns.len());
    for col in &columns {
        ctrl.extend(col.iter().map(|h| HPoint { w: 1.0, ..*h }));
    }
    NurbsSurface::from_homogeneous([p, 3], [u_knots.to_vec(), v_knots], ctrl).map_err(|e| failed(format!("closed loft: {e}")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use solidcore_record::{Vec2, Vec3};

    /// Quarter-circle profiles (two cubic pieces each) round a vertical axis, 8 of them: a
    /// hemisphere, cut into 4 faces. Its area is 2πr² (within the cubic arcs' error), and the
    /// faces meet smoothly at the seams (their normals agree there).
    #[test]
    fn closed_loft_through_radial_profiles_is_a_smooth_dome() {
        let r = 50.0;
        // A quarter circle from (r, 0) to (0, r) in two cubic pieces of 45 degrees each.
        let k = 4.0 / 3.0 * (std::f64::consts::PI / 16.0).tan();
        let piece = |a0: f64, a1: f64| {
            let (c0, s0, c1, s1) = (a0.cos(), a0.sin(), a1.cos(), a1.sin());
            Seg2::Cubic {
                p0: Vec2::new(r * c0, r * s0),
                p1: Vec2::new(r * (c0 - k * s0), r * (s0 + k * c0)),
                p2: Vec2::new(r * (c1 + k * s1), r * (s1 - k * c1)),
                p3: Vec2::new(r * c1, r * s1),
            }
        };
        let q = std::f64::consts::FRAC_PI_4;
        let profile = vec![piece(0.0, q), piece(q, 2.0 * q)];
        let sections: Vec<(Plane, Vec<Seg2>)> = (0..8)
            .map(|i| {
                let a = i as f64 * std::f64::consts::FRAC_PI_4;
                let plane = Plane { origin: Vec3::ZERO, x: Vec3::new(a.cos(), a.sin(), 0.0), y: Vec3::new(0.0, 0.0, 1.0) };
                (plane, profile.clone())
            })
            .collect();
        let dome = loft_sheet_closed(&sections, &[0, 2, 4, 6]).unwrap();
        assert_eq!(dome.face_count(), 4);
        let area = crate::measure(&dome).unwrap().area;
        let want = 2.0 * std::f64::consts::PI * r * r;
        assert!((area - want).abs() / want < 5e-3, "{area} vs {want}");
        assert!(loft_sheet_closed(&sections, &[]).is_err(), "no seams");
        assert!(loft_sheet_closed(&sections, &[9]).is_err(), "a seam past the sections");
        assert!(loft_sheet_closed(&sections[..2], &[0]).is_err(), "too few sections");
    }
}

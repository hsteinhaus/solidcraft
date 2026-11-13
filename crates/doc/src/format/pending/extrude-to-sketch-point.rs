//! An extrude's To Object (`extent.to`) can be a sketch point, `{"sketch": id, "point": "<id>"}`
//! (`p7`, `c3.center`), found in the evaluated sketch on every recompute, beside the picked
//! `[x, y, z]`. Older designs have only picks and open unchanged. The last release opens a design
//! with one by its newer format number and reports that extrude unread.

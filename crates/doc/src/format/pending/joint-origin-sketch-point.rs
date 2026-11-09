//! A joint origin can sit on a sketch point (`Snap::SketchPoint`): it saves the sketch's feature
//! id, the point's id and an optional x reference. Older designs have no such origin and need
//! nothing. The last release asks for an update instead of opening a design that has one.

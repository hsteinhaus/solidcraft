//! A component pattern is a Pattern feature that keeps the placements it made (`placements`:
//! the source occurrence, where it stood when the pattern was made, and, by place in the
//! pattern, each instance's occurrence and where the pattern last placed it), so they follow its
//! counts and spacings. Older designs have none: their patterned placements stay plain
//! placements. The last release keeps the placements and finds nothing for the feature to copy.

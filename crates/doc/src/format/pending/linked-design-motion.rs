//! A design link keeps the linked design's moving parts and joints (`DesignLink::parts`, `joints`,
//! `placements`; `LinkedOrigin::part`), so each placement moves them. Older designs have none and
//! need nothing. The last release warns about the newer format and opens the links without motion.

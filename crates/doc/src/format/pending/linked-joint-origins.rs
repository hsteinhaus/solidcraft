//! A design link keeps the joint origins its design defines (`DesignLink::origins`: name and
//! frame), and a joint origin can name one by its placement (`Snap::Linked`). Older designs have
//! neither; their links get the origins the next time they are read. The last release asks for an
//! update instead of opening a design that has them.

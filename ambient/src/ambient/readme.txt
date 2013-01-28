Ambient Caveats:

- The copy is performed by fusing two histories together (they share the same revision in one point). Therefore the element access for writing is unsafe (technically it will modify two objects at the same time).

- If revision is not reused it will consume memory until object is destroyed (rarely happends).

- The "copies" aren't cleaned up until the last object utilizing the revision will be destroyed (~1% of total objects).
  Note: If necessary it can be fixed by putting lock on a latest revision in history object. This way when one releases the lock he can safely deallocate the memory if there are no other locks. One can adapt this way for reusage by adding locks in modifying operations and removing lock in fusion but he will have to deal with double deallocation of revisions.

Ambient Caveats:

- The copy is performed by fusing two histories together (they share the same revision in one point). Therefore the element access for writing is unsafe (technically it will modify two objects at the same time).

- If revision is not reused it will consume memory until object is destroyed (rarely happends).

- The "copies" aren't cleaned up until the last object utilizing the revision will be destroyed (~1% of total objects). 

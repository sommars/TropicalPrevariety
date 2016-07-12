package Prevariety_Types is 

  type Edge;
  type Link_to_Edge is access Edge;
  type Array_of_Edge is array ( integer32 range <> ) of Link_to_Edge;

  type Edge ( ) is record
    label : integer32;
    inner_normals : ; --list of floats with some level of precision
    point_indices : ; --list of integers
    neighbor_indices : ; --list of integers
    children_indices : ; --list of integers
    parent_indices : ; --list of integers
    cone : ; --however we end up storing cones
  end record;

  package Lists_of_Edges is new Generic_Lists(Link_to_Edge);
  type Edge_List is new Lists_of_Edges.List;

--------------------------------------------------------------------------------
  type Hull;
  type Link_to_Hull is access Hull;
  type Array_of_Hull is array ( integer32 range <> ) of Link_to_Hull;

  type Hull ( ) is record
    points : ; --list of the actual points, not references to them
    index_to_point_map : ; --map the point to the index
    point_to_index_map : ; --map the index to the point
    edges : Edge_List; --list of edges sitting on the convex hull    
  end record;

  package Lists_of_Hulls is new Generic_Lists(Link_to_Hull);
  type Hull_List is new Lists_of_Hulls.List;

end Prevariety_Types;

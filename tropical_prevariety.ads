with Prevariety_Types;

package Tropical_Prevariety is

  -- Perform basic algorithm where all cones are intersected.
  function Common_Refinement(Hulls: Prevariety_Types.Hull_List) 
    return Prevariety_Types.Cone_List;

end Tropical_Prevariety;

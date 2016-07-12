with Prevariety_Types;

package body Tropical_Prevariety is 

  -- Perform basic algorithm where all cones are intersected.
  function Common_Refinement(Hulls: Prevariety_Types.Hull_List) 
  return Prevariety_Types.Cone_List is

    cs, new_cs : Prevariety_Types.Cone_Set;

  begin
    --TODO: somehow initialize the cone set
    
    --Possibly indexing into Hulls incorrectly.
    for i in Hulls'length loop
      if i = 1 then
        for e in Hulls'first'edges loop
          --TODO: somehow add e'cone to each set
        end loop;
      else
        for e in Hulls(i)'edges loop
          --TODO: somehow intersect e'cone with each element in cs. add to new_cs
        end loop;
        cs := new_cs;
      end if;
    end loop;
    return cs;
  end Common_Refinement;

end Tropical_Prevariety;

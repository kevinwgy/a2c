#ifndef _MATERIAL_OPERATOR_H_
#define _MATERIAL_OPERATOR_H_
#include<IoData.h>
#include<string>
#include<vector>
#include<map>
#include<cassert>

/***************************************************************
 * class MaterialOperator is responsible for performing material-
 * dependent calculations
 **************************************************************/

class MaterialOperator
{

  int material_id; //!< my ID
  MaterialData *iod_mat_ptr; //!< input variables for THIS MATERIAL

  std::vector<std::pair<int, std::string> > species_l2g; //!< <global ID, name> of my species
  std::map<int, int> species_g2l; //!< global ID --> local species index

  double xmin;

public:

  MaterialOperator();
  ~MaterialOperator();

  void Setup(int material_id_, MaterialData &iod_mat, std::vector<string> &global_species);

  inline int GetMaterialID() {return material_id;}
  inline int GetNumberOfSpecies() {return species_l2g.size();}

  inline double GetMinMolarFraction() {return xmin;}

  inline int LocalToGlobalSpeciesID(int local_id)
    {assert(local_id>=0 && local_id<(int)species_l2g.size()); return species_l2g[local_id].first;}
  inline int GlobalToLocalSpeciesID(int global_id)
    {auto it = species_g2l.find(global_id); return it==species_g2l.end() ? -1 : it->second;}
  inline string GetSpeciesNameByLocalID(int local_id)
    {assert(local_id>=0 && local_id<(int)species_l2g.size()); return species_l2g[local_id].second;}

};

#endif

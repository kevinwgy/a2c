#include<MaterialOperator.h>
#include<Utils.h>
#include<set>

//------------------------------------------------------------------------

MaterialOperator::MaterialOperator() 
                : material_id(-1), iod_mat_ptr(NULL), xmin(-1.0)
{ }

//------------------------------------------------------------------------

MaterialOperator::~MaterialOperator()
{ }

//------------------------------------------------------------------------

void
MaterialOperator::Setup(int material_id_, MaterialData &iod_mat, vector<string> &global_species)
{
  material_id = material_id_;
  iod_mat_ptr = &iod_mat;

  species_l2g.resize(iod_mat_ptr->nSpecies, std::make_pair(-1,""));
  std::set<int> checker;
  for(int i=0; i<iod_mat_ptr->nSpecies; i++) {
    int gid = iod_mat_ptr->species[i];
    if(gid<0 || gid>=(int)global_species.size()) {
      print_error("*** Error: Material[%d] has undefined species (id=%d).\n",
                  material_id, gid);
      exit_mpi();
    }
    species_l2g[i]   = std::make_pair(gid, global_species[gid]);
    species_g2l[gid] = i;
    checker.insert(gid);
  }
  if((int)checker.size() != iod_mat_ptr->nSpecies) {
    print_error("*** Error: Detected duplicate species IDs in Material[%d].\n", material_id); 
    exit_mpi();
  }

  xmin = iod_mat_ptr->min_molar_fraction;
  if(xmin<0.0) {
    print_error("*** Error: Detected negative molar fraction (%e) in Material[%d].\n",
                xmin, material_id);
    exit_mpi();
  }

}

//------------------------------------------------------------------------



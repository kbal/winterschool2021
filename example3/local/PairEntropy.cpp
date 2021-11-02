/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2017 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_tools_NeighborListParallel_h
#define __PLUMED_tools_NeighborListParallel_h

#include "tools/Vector.h"
#include "tools/AtomNumber.h"
#include "tools/Communicator.h"
#include "tools/Log.h"

#include <vector>

namespace PLMD {

class Pbc;

/// \ingroup TOOLBOX
/// A class that implements neighbor lists from two lists or a single list of atoms
class NeighborListParallel
{
  bool do_pair_,do_pbc_,twolists_,firsttime_;
/// Choose between half and full list
  bool do_full_list_;
  const PLMD::Pbc* pbc_;
/// Full list of atoms involved
  std::vector<PLMD::AtomNumber> fullatomlist_,requestlist_;
/// Positions during the building of the NL
/// They are used for checking if the next build is dangerous
  std::vector<PLMD::Vector> positions_old_;
/// Neighbor list
  std::vector<std::vector<unsigned>> neighbors_;
  //std::vector<std::pair<unsigned,unsigned> > neighbors_;
/// List of the "local" atoms each thread has
  std::vector<unsigned> local_atoms_;
  double distance_, skin_;
  int stride_;
  unsigned nlist0_,nlist1_,nallpairs_,lastupdate_;
  unsigned dangerousBuilds_, numberOfBuilds_;
  double maxLoadImbalance_,avgLoadImbalance_, avgTotalNeighbors_;
/// Return the pair of indexes in the positions array
/// of the two atoms forming the i-th pair among all possible pairs
  //std::pair<unsigned,unsigned> getIndexPair(unsigned i);
/// Communicator
  Communicator& mycomm;
/// Log
  Log& mylog;
/// MPI stuff
  unsigned mpi_rank, mpi_stride;
/// Linked list stuff
  double binsizex, binsizey, binsizez, bininvx, bininvy, bininvz ;
  bool do_linked_list_;
/// Number of bins in x,y,z
  int nbinx, nbiny, nbinz;
/// Head of chain for each bin
  std::vector<int> binhead;
/// Linked list for each bin
  std::vector<unsigned> bins;
/// sx,sy,sz = max range of stencil in each dim
  int sx,sy,sz;
/// Number of stencils
  unsigned nstencil;
/// Stencil indices
  std::vector<unsigned> stencilx, stencily, stencilz;
public:
  NeighborListParallel(const std::vector<PLMD::AtomNumber>& list0,
               const std::vector<PLMD::AtomNumber>& list1,
               const bool& do_pair, const bool& do_pbc, const PLMD::Pbc& pbc, Communicator& cc,
               Log& log, const double& distance=1.0e+30,  const bool& do_full_list=false, const int& stride=0, const double& skin=0.1);
  NeighborListParallel(const std::vector<PLMD::AtomNumber>& list0, const bool& do_pbc,
               const PLMD::Pbc& pbc, Communicator& cc, Log& log, const double& distance=1.0e+30,
               const bool& do_full_list=false, const int& stride=0, const double& skin=0.1);
/// Return the list of all atoms. These are needed to rebuild the neighbor list.
  std::vector<PLMD::AtomNumber>& getFullAtomList();
/// Check if the nieghbor list must be rebuilt
  bool isListStillGood(const std::vector<Vector>& positions);
/// Update the indexes in the neighbor list to match the
/// ordering in the new positions array
/// and return the new list of atoms that must be requested to the main code
//  std::vector<PLMD::AtomNumber>& getReducedAtomList();
/// Update the neighbor list and prepare the new
/// list of atoms that will be requested to the main code
  void update(const std::vector<PLMD::Vector>& positions);
/// Construct a half neighbor list, containing only i,j and not j,i
  void updateHalfList(const std::vector<Vector>& positions);
/// Construct a full neighbor list, containing both i,j and j,i
  void updateFullList(const std::vector<Vector>& positions);
/// Construct a full neighbor list, containing both i,j and j,i ; Linked lists are used
  void updateFullListWithLinkedList(const std::vector<Vector>& positions);
/// Construct a half neighbor list, containing only i,j and not j,i ; Linked lists are used
  void updateHalfListWithLinkedList(const std::vector<Vector>& positions);
/// Get the update stride of the neighbor list
  int getStride() const;
/// Get the last step in which the neighbor list was updated
  unsigned getLastUpdate() const;
/// Set the step of the last update
  void setLastUpdate(unsigned step);
/// Get the size of the neighbor list - Obsolete
  unsigned size() const;
/// Get the i-th pair of the neighbor list - Obsolete
  std::pair<unsigned,unsigned> getClosePair(unsigned i) const;
/// Get the list of neighbors of the i-th atom
  std::vector<unsigned> getNeighbors(unsigned i);
  ~NeighborListParallel() {}
/// Print statistics of neighbor list
  void printStats();
/// Gather statistics of neighbor list
  void gatherStats(const std::vector<PLMD::Vector>& positions);
/// Get the local number of atoms (the ones each thread is handling)
  unsigned getNumberOfLocalAtoms() const;
/// Get the index if the i-th local atom
  unsigned getIndexOfLocalAtom(unsigned i) const;
/// Prepare linked lists
  void prepareLinkedList();
/// Coordinate to bin
  unsigned coord2bin(const Vector& position);
/// Coordinate to bin, get also position in x,y,z
  unsigned coord2bin(const Vector& position, unsigned& ibinx, unsigned& ibiny, unsigned& ibinz);
/// Find the k-th neighboring cell of cell atombinx,atombiny,atombinz
  unsigned neighborCell(const unsigned& k, const unsigned& atombinx, const unsigned& atombiny, const unsigned& atombinz);
/// Calculate the distance with or without pbc, as defined with the do_pbc_ flag
  Vector distance(const Vector& position1, const Vector& position2);
};

}

#endif
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2017 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "NeighborListParallel.h"
#include "tools/Vector.h"
#include "tools/Pbc.h"
#include "tools/AtomNumber.h"
#include "tools/Tools.h"
#include <vector>
#include <algorithm>
#include "tools/Communicator.h"
#include "tools/Log.h"

namespace PLMD {
using namespace std;

NeighborListParallel::NeighborListParallel(const vector<AtomNumber>& list0, const vector<AtomNumber>& list1,
                           const bool& do_pair, const bool& do_pbc, const Pbc& pbc, Communicator& cc, Log& log,
                           const double& distance, const bool& do_full_list, const int& stride, const double& skin): 
  do_pair_(do_pair), do_pbc_(do_pbc), pbc_(&pbc),
  distance_(distance), stride_(stride), mycomm(cc), mylog(log),
  skin_(skin), do_full_list_(do_full_list)
{
// store full list of atoms needed
  fullatomlist_=list0;
  fullatomlist_.insert(fullatomlist_.end(),list1.begin(),list1.end());
  nlist0_=list0.size();
  nlist1_=list1.size();
  twolists_=true;
  if(!do_pair) {
    nallpairs_=nlist0_*nlist1_;
  } else {
    plumed_assert(nlist0_==nlist1_);
    nallpairs_=nlist0_;
  }
  lastupdate_=0;
  positions_old_.resize(fullatomlist_.size());
  dangerousBuilds_=0;
  firsttime_=true;
  numberOfBuilds_=0;
  avgTotalNeighbors_=0.;
  maxLoadImbalance_=2.;
  avgLoadImbalance_=0.;
  // Setup mpi
  mpi_rank=mycomm.Get_rank();
  mpi_stride=mycomm.Get_size();
}

NeighborListParallel::NeighborListParallel(const vector<AtomNumber>& list0, const bool& do_pbc,
                           const Pbc& pbc, Communicator& cc, Log& log, const double& distance,
                           const bool& do_full_list, const int& stride, const double& skin):
  do_pbc_(do_pbc), pbc_(&pbc),
  distance_(distance), stride_(stride), mycomm(cc), mylog(log),
  skin_(skin) , do_full_list_(do_full_list)
{
  fullatomlist_=list0;
  nlist0_=list0.size();
  twolists_=false;
  nallpairs_=nlist0_*(nlist0_-1)/2;
  lastupdate_=0;
  positions_old_.resize(fullatomlist_.size());
  dangerousBuilds_=0;
  firsttime_=true;
  numberOfBuilds_=0;
  avgTotalNeighbors_=0.;
  maxLoadImbalance_=2.;
  avgLoadImbalance_=0.;
  // Setup mpi
  mpi_rank=mycomm.Get_rank();
  mpi_stride=mycomm.Get_size();
}

vector<AtomNumber>& NeighborListParallel::getFullAtomList() {
  return fullatomlist_;
}

bool NeighborListParallel::isListStillGood(const vector<Vector>& positions) {
  bool flag=true;
  plumed_assert(positions.size()==fullatomlist_.size());
  for(unsigned int i=0;i<fullatomlist_.size();i++) {
    Vector distance;
    if(do_pbc_) {
       distance=pbc_->distance(positions[i],positions_old_[i]);
    } else {
       distance=delta(positions[i],positions_old_[i]);
    }
    if (modulo(distance)>skin_) {
       flag=false;
       break;
    }
  }
  return flag;
}

void NeighborListParallel::printStats() {
  mylog.printf("Neighbor list statistics\n");
  mylog.printf("Total # of neighbors = %f \n", avgTotalNeighbors_);
  mylog.printf("Ave neighs/atom = %f \n", avgTotalNeighbors_ /(double) nlist0_);
  mylog.printf("Neighbor list builds = %d \n",numberOfBuilds_);
  mylog.printf("Dangerous builds = %d \n",dangerousBuilds_);
  mylog.printf("Average load imbalance (min/max) = %f \n",avgLoadImbalance_);
  mylog.printf("Maximum load imbalance (min/max) = %f \n",maxLoadImbalance_);
  if (do_linked_list_) {
     mylog.printf("Number of bins in linked list = %d %d %d \n", nbinx, nbiny, nbinz);
  }
}

void NeighborListParallel::update(const vector<Vector>& positions) {
  // clear previous list
  neighbors_.clear();
  local_atoms_.clear();
  // check if positions array has the correct length
  plumed_assert(positions.size()==fullatomlist_.size());
  // Prepare linked lists
  prepareLinkedList();
  // Decide whether to do linked lists or the N^2 calculation
  do_linked_list_=true;
  if ((2*sx+1) >= nbinx || (2*sy+1) >= nbiny || (2*sz+1) >= nbinz) do_linked_list_=false;
  if (do_linked_list_) {
     if (!do_full_list_) updateHalfListWithLinkedList(positions);
     else updateFullListWithLinkedList(positions);
  } else {
     if (!do_full_list_) updateHalfList(positions);
     else updateFullList(positions);
  }
  gatherStats(positions);
  // Store positions for checking
  for(unsigned int i=0;i<fullatomlist_.size();i++) {
     positions_old_[i]=positions[i];
  }
}

void NeighborListParallel::updateFullListWithLinkedList(const vector<Vector>& positions) {
  const double d2=distance_*distance_;
  // Set bin head of chain to -1.0
  binhead = std::vector<int>(nbinx*nbiny*nbinz,-1.0);
  if (!twolists_) {
    neighbors_.resize(nlist0_);
    // Construct linked list and assign head of chain
    bins.resize(nlist0_);
    for(unsigned int i=0;i<nlist0_;i+=1) {
       unsigned ibin = coord2bin(positions[i]);
       bins[i] = binhead[ibin];
       binhead[ibin] = i;
    }
    // Calculate
    for(unsigned int i=mpi_rank;i<nlist0_;i+=mpi_stride) {
       local_atoms_.push_back(i);
       Vector position_i=positions[i];
       unsigned atombinx, atombiny, atombinz, atombin;
       atombin = coord2bin(positions[i],atombinx, atombiny, atombinz);
       for(unsigned k=0; k < nstencil; k++) {
          // Loop over atoms in the k-th neighboring cell
          unsigned kbin = neighborCell(k,atombinx, atombiny, atombinz) ;
          for (int j = binhead[kbin]; j >= 0; j = bins[j]) {
             double value=modulo2(distance(position_i,positions[j]));
             if(value<=d2) neighbors_[i].push_back(j);
          }
       }
    }
  } else if (twolists_ && !do_pair_) {
    neighbors_.resize(nlist0_+nlist1_);
    // Construct linked list and assign head of chain
    bins.resize(nlist1_);
    for(unsigned int i=0;i<nlist1_;i++) {
       unsigned ibin = coord2bin(positions[i+nlist0_]);
       bins[i] = binhead[ibin];
       binhead[ibin] = i;
    }
    // Calculate
    for(unsigned int i=mpi_rank;i<nlist0_;i+=mpi_stride) {
       local_atoms_.push_back(i);
       Vector position_i=positions[i];
       unsigned atombinx, atombiny, atombinz, atombin;
       atombin = coord2bin(positions[i],atombinx, atombiny, atombinz);
       for(unsigned k=0; k < nstencil; k++) {
          // Loop over atoms in the k-th neighboring cell
          unsigned kbin = neighborCell(k,atombinx, atombiny, atombinz) ;
          for (int j = binhead[kbin]; j >= 0; j = bins[j]) {
             double value=modulo2(distance(position_i,positions[j+nlist0_]));
             if(value<=d2) {
                neighbors_[i].push_back(j+nlist0_);
                neighbors_[j+nlist0_].push_back(i);
             }
          }
       }
    }
  }
}

void NeighborListParallel::updateHalfListWithLinkedList(const vector<Vector>& positions) {
  const double d2=distance_*distance_;
  neighbors_.resize(nlist0_);
  // Set bin head of chain to -1.0
  binhead = std::vector<int>(nbinx*nbiny*nbinz,-1.0);
  if (!twolists_) {
    // Construct linked list and assign head of chain
    bins.resize(nlist0_);
    for(unsigned int i=0;i<nlist0_;i+=1) {
       unsigned ibin = coord2bin(positions[i]);
       bins[i] = binhead[ibin];
       binhead[ibin] = i;
    }
    // Calculate
    for(unsigned int i=mpi_rank;i<nlist0_;i+=mpi_stride) {
       local_atoms_.push_back(i);
       Vector position_i=positions[i];
       unsigned atombinx, atombiny, atombinz, atombin;
       atombin = coord2bin(positions[i],atombinx, atombiny, atombinz);
       for(unsigned k=0; k < nstencil; k++) {
          // Loop over atoms in the k-th neighboring cell
          unsigned kbin = neighborCell(k,atombinx, atombiny, atombinz) ;
          for (int j = binhead[kbin]; j >= 0; j = bins[j]) {
             if (j>=i) continue;
             double value=modulo2(distance(position_i,positions[j]));
             if(value<=d2) neighbors_[i].push_back(j);
          }
       }
    }
  } else if (twolists_ && !do_pair_) {
    // Construct linked list and assign head of chain
    bins.resize(nlist1_);
    for(unsigned int i=0;i<nlist1_;i++) {
       unsigned ibin = coord2bin(positions[i+nlist0_]);
       bins[i] = binhead[ibin];
       binhead[ibin] = i;
    }
    // Calculate
    for(unsigned int i=mpi_rank;i<nlist0_;i+=mpi_stride) {
       local_atoms_.push_back(i);
       Vector position_i=positions[i];
       unsigned atombinx, atombiny, atombinz, atombin;
       atombin = coord2bin(positions[i],atombinx, atombiny, atombinz);
       for(unsigned k=0; k < nstencil; k++) {
          // Loop over atoms in the k-th neighboring cell
          unsigned kbin = neighborCell(k,atombinx, atombiny, atombinz) ;
          for (int j = binhead[kbin]; j >= 0; j = bins[j]) {
             double value=modulo2(distance(position_i,positions[j+nlist0_]));
             if(value<=d2) {
                neighbors_[i].push_back(j+nlist0_);
             }
          }
       }
    }
  }
}


void NeighborListParallel::updateFullList(const vector<Vector>& positions) {
  const double d2=distance_*distance_;
  if (!twolists_) {
    neighbors_.resize(nlist0_);
    for(unsigned int i=mpi_rank;i<nlist0_;i+=mpi_stride) {
       local_atoms_.push_back(i);
       Vector position_i=positions[i];
       for(unsigned int j=0;j<nlist0_;j+=1) {
          double value=modulo2(distance(position_i,positions[j]));
          if(value<=d2) neighbors_[i].push_back(j);
       }
    }
  } else if(twolists_ && do_pair_) {
    neighbors_.resize(nlist0_+nlist1_);
    for(unsigned int i=mpi_rank;i<nlist0_;i+=mpi_stride) {
       local_atoms_.push_back(i);
       double value=modulo2(distance(positions[i],positions[nlist0_+i]));
       if(value<=d2) {
          neighbors_[i].push_back(nlist0_+i);
          neighbors_[nlist0_].push_back(i);
       }
    }
  } else if (twolists_ && !do_pair_) {
    neighbors_.resize(nlist0_+nlist1_);
    for(unsigned int i=mpi_rank;i<nlist0_;i+=mpi_stride) {
       local_atoms_.push_back(i);
       Vector position_i=positions[i];
       for(unsigned int j=0;j<nlist1_;j+=1) {
          double value=modulo2(distance(position_i,positions[nlist0_+j]));
          if(value<=d2) {
             neighbors_[i].push_back(nlist0_+j);
             neighbors_[nlist0_+j].push_back(i);
          }
       }
    }
  }
}

void NeighborListParallel::updateHalfList(const vector<Vector>& positions) {
  const double d2=distance_*distance_;
  neighbors_.resize(nlist0_);
  if (!twolists_) {
    for(unsigned int i=mpi_rank;i<nlist0_;i+=mpi_stride) {
       local_atoms_.push_back(i);
       Vector position_i=positions[i];
       for(unsigned int j=i+1;j<nlist0_;j+=1) {
          double value=modulo2(distance(position_i,positions[j]));
          if(value<=d2) neighbors_[i].push_back(j);
       }
    }
  } else if(twolists_ && do_pair_) {
    for(unsigned int i=mpi_rank;i<nlist0_;i+=mpi_stride) {
       local_atoms_.push_back(i);
       double value=modulo2(distance(positions[i],positions[nlist0_+i]));
       if(value<=d2) neighbors_[i].push_back(nlist0_+i);
    }
  } else if (twolists_ && !do_pair_) {
    for(unsigned int i=mpi_rank;i<nlist0_;i+=mpi_stride) {
       local_atoms_.push_back(i);
       Vector position_i=positions[i];
       for(unsigned int j=0;j<nlist1_;j++) {
          double value=modulo2(distance(position_i,positions[nlist0_+j]));
          if(value<=d2) {
             neighbors_[i].push_back(nlist0_+j);
          }
       }
    }
  }
}

void NeighborListParallel::gatherStats(const vector<Vector>& positions) {
  // Check if rebuilt was dangerous
  if (!firsttime_ && !isListStillGood(positions)) {
     dangerousBuilds_++;
  }
  firsttime_=false;
  numberOfBuilds_++;
  unsigned neighNum=0;
  for(unsigned int i=0;i<nlist0_;i++) {
     neighNum += neighbors_[i].size();
  }
  unsigned allNeighNum=0;
  std::vector<unsigned> neighbors_ranks_(mycomm.Get_size());
  mycomm.Allgather(&neighNum,1,&neighbors_ranks_[0],1);
  for(unsigned int i=0;i<mycomm.Get_size();i+=1) allNeighNum+=neighbors_ranks_[i];
  auto min_element_ = *std::min_element(neighbors_ranks_.begin(), neighbors_ranks_.end());
  auto max_element_ = *std::max_element(neighbors_ranks_.begin(), neighbors_ranks_.end());
  double loadImbalance=min_element_ / (double) max_element_;
  if (maxLoadImbalance_>loadImbalance) maxLoadImbalance_=loadImbalance;
  avgLoadImbalance_ += (loadImbalance-avgLoadImbalance_)/numberOfBuilds_;
  avgTotalNeighbors_ += (allNeighNum-avgTotalNeighbors_)/numberOfBuilds_;
}

unsigned NeighborListParallel::getNumberOfLocalAtoms() const {
  return local_atoms_.size();
}

unsigned NeighborListParallel::getIndexOfLocalAtom(unsigned i) const {
  return local_atoms_[i];
}

int NeighborListParallel::getStride() const {
  return stride_;
}

unsigned NeighborListParallel::getLastUpdate() const {
  return lastupdate_;
}

void NeighborListParallel::setLastUpdate(unsigned step) {
  lastupdate_=step;
}

// This function is obsolete and should be changed
unsigned NeighborListParallel::size() const {
  return neighbors_.size();
}

// This function is obsolete and should be changed
pair<unsigned,unsigned> NeighborListParallel::getClosePair(unsigned i) const {
  return pair<unsigned,unsigned>(i,i);
}

vector<unsigned> NeighborListParallel::getNeighbors(unsigned index) {
  return neighbors_[index];
}

Vector NeighborListParallel::distance(const Vector& position1, const Vector& position2) {
  Vector distance;
  if(do_pbc_) {
    distance=pbc_->distance(position1,position2);
  } else {
    distance=delta(position1,position2);
  }
  return distance;
}

// Linked lists stuff

void NeighborListParallel::prepareLinkedList() {
  // Determine optimal number of cells in each direction
  double binsize_optimal = 0.5*distance_;
  double binsizeinv = 1.0/binsize_optimal;
  Tensor bbox=pbc_->getBox();
  nbinx = static_cast<int> (bbox[0][0]*binsizeinv);
  nbiny = static_cast<int> (bbox[1][1]*binsizeinv);
  nbinz = static_cast<int> (bbox[2][2]*binsizeinv);
  if (nbinx == 0) nbinx = 1;
  if (nbiny == 0) nbiny = 1;
  if (nbinz == 0) nbinz = 1;
  // Determine cell sizes
  binsizex = 1.0 / nbinx;
  binsizey = 1.0 / nbiny;
  binsizez = 1.0 / nbinz;
  bininvx = 1.0 / binsizex;
  bininvy = 1.0 / binsizey;
  bininvz = 1.0 / binsizez;
  // sx,sy,sz = max range of stencil in each dim
  Vector scaled_distance_=pbc_->realToScaled(Vector(distance_,distance_,distance_));
  sx = static_cast<int> (scaled_distance_[0]*bininvx);
  if (sx*binsizex < scaled_distance_[0]) sx++;
  sy = static_cast<int> (scaled_distance_[1]*bininvy);
  if (sy*binsizey < scaled_distance_[1]) sy++;
  sz = static_cast<int> (scaled_distance_[2]*bininvz);
  if (sz*binsizez < scaled_distance_[2]) sz++;
  // Create stencil
  nstencil = (2*sx+1)*(2*sy+1)*(2*sz+1);
  stencilx.resize(nstencil);
  stencily.resize(nstencil);
  stencilz.resize(nstencil);
  unsigned k=0;
  for(int ix=-sx;ix<=sx;ix++) {
     for(int iy=-sy;iy<=sy;iy++) {
        for(int iz=-sz;iz<=sz;iz++) {
           stencilx[k]=ix;
           stencily[k]=iy;
           stencilz[k]=iz;
           k++;
        }
     }
  }
   // clear previous list
  bins.clear();
  binhead.clear();
}

unsigned NeighborListParallel::neighborCell(const unsigned& k, const unsigned& atombinx, const unsigned& atombiny, const unsigned& atombinz) {
  int ibinx = atombinx+stencilx[k];
  if (ibinx<0) ibinx += nbinx; 
  if (ibinx>=nbinx) ibinx -= nbinx;
  int ibiny = atombiny+stencily[k];
  if (ibiny<0) ibiny += nbiny; 
  if (ibiny>=nbiny) ibiny -= nbiny;
  int ibinz = atombinz+stencilz[k];
  if (ibinz<0) ibinz += nbinz; 
  if (ibinz>=nbinz) ibinz -= nbinz;
  unsigned neighbor = ibinx + ibiny*nbinx + ibinz*nbinx*nbiny;
  return neighbor; 
}

unsigned NeighborListParallel::coord2bin(const Vector& position) {
  Vector pbc_position=pbc_->distance(Vector(0.0,0.0,0.0),position);
  Vector scaled_position=pbc_->realToScaled(pbc_position);
  if (scaled_position[0]<0.) scaled_position[0] += 1;
  if (scaled_position[1]<0.) scaled_position[1] += 1;
  if (scaled_position[2]<0.) scaled_position[2] += 1;
  unsigned ibinx, ibiny, ibinz, ibin;
  ibinx = static_cast<int> (scaled_position[0]*bininvx);
  ibiny = static_cast<int> (scaled_position[1]*bininvy);
  ibinz = static_cast<int> (scaled_position[2]*bininvz);
  ibin = ibinx + ibiny*nbinx + ibinz*nbinx*nbiny;
  return ibin;
}

unsigned NeighborListParallel::coord2bin(const Vector& position, unsigned& ibinx, unsigned& ibiny, unsigned& ibinz) {
  Vector pbc_position=pbc_->distance(Vector(0.0,0.0,0.0),position);
  Vector scaled_position=pbc_->realToScaled(pbc_position);
  if (scaled_position[0]<0.) scaled_position[0] += 1;
  if (scaled_position[1]<0.) scaled_position[1] += 1;
  if (scaled_position[2]<0.) scaled_position[2] += 1;
  unsigned ibin;
  ibinx = static_cast<int> (scaled_position[0]*bininvx);
  ibiny = static_cast<int> (scaled_position[1]*bininvy);
  ibinz = static_cast<int> (scaled_position[2]*bininvz);
  ibin = ibinx + ibiny*nbinx + ibinz*nbinx*nbiny;
  return ibin;
}

}
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2016 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "colvar/Colvar.h"
#include "colvar/ActionRegister.h"
#include "NeighborListParallel.h"
#include "tools/Communicator.h"
#include "tools/Tools.h"
#include "tools/IFile.h"
#include "math.h"

#include <string>

using namespace std;

namespace PLMD{
namespace colvar{

//+PLUMEDOC COLVAR PAIRENTROPY
/*
Calculate the global pair entropy using the expression:
\f[
s=-2\pi\rho k_B \int\limits_0^{r_{\mathrm{max}}} \left [ g(r) \ln g(r) - g(r) + 1 \right ] r^2 dr .
\f]
where \f$ g(r) $\f is the pair distribution function and \f$ r_{\mathrm{max}} $\f is a cutoff in the integration (MAXR).
For the integration the interval from 0 to  \f$ r_{\mathrm{max}} $\f is partitioned in NHIST equal intervals. 
To make the calculation of \f$ g(r) $\f differentiable, the following function is used:
\f[
g(r) = \frac{1}{4 \pi \rho r^2} \sum\limits_{j} \frac{1}{\sqrt{2 \pi \sigma^2}} e^{-(r-r_{ij})^2/(2\sigma^2)} ,
\f]
where \f$ \rho $\f is the density and \f$ sigma $\f is a broadening parameter (SIGMA).  
\par Example)
The following input tells plumed to calculate the pair entropy of atoms 1-250 with themselves.
\verbatim
PAIRENTROPY ...
 LABEL=s2
 ATOMS=1-250
 MAXR=0.65
 SIGMA=0.01
 NLIST
 NL_CUTOFF=0.75
 NL_STRIDE=10
... PAIRENTROPY
\endverbatim
*/
//+ENDPLUMEDOC

class PairEntropy : public Colvar {
  bool pbc, serial;
  // Neighbor list stuff
  bool doneigh;
  NeighborListParallel *nl;
  vector<AtomNumber> atoms_lista;
  bool invalidateList;
  bool firsttime;
  // Others
  double maxr, sigma;
  unsigned nhist;
  double rcut2;
  double sqrt2piSigma, sigmaSqr2, sigmaSqr;
  double deltar;
  unsigned deltaBin;
  double density_given;
  std::vector<double> vectorX, vectorX2;
  // Integration routines
  double integrate(vector<double> integrand, double delta)const;
  Vector integrate(vector<Vector> integrand, double delta)const;
  Tensor integrate(vector<Tensor> integrand, double delta)const;
  // Kernel to calculate g(r)
  double kernel(double distance, double invNormKernel, double&der)const;
  // Output gofr and integrand
  bool doOutputGofr;
  bool doOutputIntegrand;
  unsigned outputStride;
  void outputGofr(vector<double> gofr);
  void outputIntegrand(vector<double> integrand);
  mutable PLMD::OFile gofrOfile, integrandOfile;
  // Reference g(r)
  bool doReferenceGofr;
  bool noPrefactor;
  std::vector<double> referenceGofr;
  // Average g(r)
  bool doAverageGofr;
  vector<double> avgGofr;
  unsigned iteration;
  // Low communication variant
  bool doLowComm;
public:
  explicit PairEntropy(const ActionOptions&);
  ~PairEntropy();
  virtual void calculate();
  virtual void prepare();
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(PairEntropy,"PAIRENTROPY")

void PairEntropy::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords(keys);
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.addFlag("NLIST",false,"Use a neighbour list to speed up the calculation");
  keys.addFlag("OUTPUT_GOFR",false,"Output g(r)");
  keys.addFlag("OUTPUT_INTEGRAND",false,"Output integrand");
  keys.add("optional","OUTPUT_STRIDE","The frequency with which the output is written to files");
  keys.addFlag("AVERAGE_GOFR",false,"Average g(r) over time");
  keys.addFlag("NO_PREFACTOR",false,"Do not include the prefactor -4*pi*density");
  keys.add("optional","NL_CUTOFF","The cutoff for the neighbor list");
  keys.add("optional","NL_STRIDE","The frequency with which we are updating the atoms in the neighbour list. If non specified or negative, it checks every step and rebuilds as needed.");
  keys.add("optional","DENSITY","Density to normalize the g(r). If not specified, N/V is used");
  keys.add("atoms","ATOMS","List of atoms");
  keys.add("compulsory","MAXR","1.","Maximum distance for the radial distribution function ");
  keys.add("optional","NHIST","Number of bins in the rdf ");
  keys.add("compulsory","SIGMA","0.01","Width of gaussians ");
  keys.add("optional","REFERENCE_GOFR_FNAME","the name of the file with the reference g(r)");
  keys.addFlag("LOW_COMM",false,"Use an algorithm with less communication between processors");
}

PairEntropy::PairEntropy(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
pbc(true),
serial(false),
invalidateList(true),
firsttime(true)
{

  parseFlag("SERIAL",serial);

  parseAtomList("ATOMS",atoms_lista);

  log.printf("  using periodic boundary conditions\n");

// neighbor list stuff
  doneigh=false;
  bool nl_full_list=false;
  double nl_cut=0.0;
  double nl_skin;
  int nl_st=-1;
  parseFlag("NLIST",doneigh);
  if(doneigh){
   parse("NL_CUTOFF",nl_cut);
   if(nl_cut<=0.0) error("NL_CUTOFF should be explicitly specified and positive");
   parse("NL_STRIDE",nl_st);
  }

  density_given = -1;
  parse("DENSITY",density_given);
  if (density_given>0) log.printf("  The g(r) will be normalized with a density %f . \n", density_given);
  else log.printf("  The g(r) will be normalized with a density N/V . \n");

  addValueWithDerivatives(); setNotPeriodic();

  parse("MAXR",maxr);
  log.printf("  Integration in the interval from 0. to %f \n", maxr );
  parse("SIGMA",sigma);
  log.printf("  The pair distribution function is calculated with a Gaussian kernel with deviation %f \n", sigma);
  double rcut = maxr + 3*sigma;
  rcut2 = (maxr + 3*sigma)*(maxr + 3*sigma);  // 3*sigma is hard coded
  if(doneigh){
    if(nl_cut<rcut) error("NL_CUTOFF should be larger than MAXR + 3*SIGMA");
    nl_skin=nl_cut-(maxr+2*sigma);
  }
  nhist=ceil(maxr/sigma) + 1; // Default value
  parse("NHIST",nhist);
  log.printf("  The interval is partitioned in %u equal parts and the integration is perfromed with the trapezoid rule. \n", nhist );
 
  doOutputGofr=false;
  parseFlag("OUTPUT_GOFR",doOutputGofr);
  if (doOutputGofr) { 
     log.printf("  The g(r) will be written to a file \n.");
     gofrOfile.link(*this);
     gofrOfile.open("gofr.txt");
  }
  doOutputIntegrand=false;
  parseFlag("OUTPUT_INTEGRAND",doOutputIntegrand);
  if (doOutputIntegrand) {
     log.printf("  The integrand will be written to a file \n.");
     integrandOfile.link(*this);
     integrandOfile.open("integrand.txt");
  }
  outputStride=1;
  parse("OUTPUT_STRIDE",outputStride);
  if (outputStride!=1 && !doOutputGofr && !doOutputIntegrand) error("Cannot specify OUTPUT_STRIDE if OUTPUT_GOFR or OUTPUT_INTEGRAND not used");
  if (outputStride<1) error("The output stride specified with OUTPUT_STRIDE must be greater than or equal to one.");
  if (outputStride>1) log.printf("  The output stride to write g(r) or the integrand is %d \n", outputStride);

  doReferenceGofr=false;
  std::string referenceGofrFileName;
  parse("REFERENCE_GOFR_FNAME",referenceGofrFileName); 
  if (!referenceGofrFileName.empty() ) {
    log.printf("  Reading a reference g(r) from the file %s . \n", referenceGofrFileName.c_str() );
    doReferenceGofr=true;
    IFile ifile; 
    ifile.link(*this);
    ifile.open(referenceGofrFileName);
    referenceGofr.resize(nhist);
    for(unsigned int i=0;i<nhist;i++) {
       double tmp_r;
       ifile.scanField("r",tmp_r).scanField("gofr",referenceGofr[i]).scanField();
    }
  }

  noPrefactor=false;
  parseFlag("NO_PREFACTOR",noPrefactor);

  doAverageGofr=false;
  parseFlag("AVERAGE_GOFR",doAverageGofr);
  if (doAverageGofr) {
     iteration = 1;
     log.printf("  The g(r) will be averaged over all frames");
     avgGofr.resize(nhist);
  }

  doLowComm=false;
  parseFlag("LOW_COMM",doLowComm);
  if (doLowComm) {
     log.printf("  Using the low communication variant of the algorithm");
     nl_full_list=true;
     if (!doneigh) error("LOW_COMM can only be used with neighbor lists");
  }

  checkRead();

  // Setup neighbor list
  if (doneigh) {
    nl= new NeighborListParallel(atoms_lista,pbc,getPbc(),comm,log,nl_cut,nl_full_list,nl_st,nl_skin);
    requestAtoms(nl->getFullAtomList());
    log.printf("  using neighbor lists with\n");
    log.printf("  cutoff %f, and skin %f\n",nl_cut,nl_skin);
    if(nl_st>=0){
      log.printf("  update every %d steps\n",nl_st);
    } else {
      log.printf("  checking every step for dangerous builds and rebuilding as needed\n");
    }
    if (nl_full_list) {
      log.printf("  using a full neighbor list\n");
    } else {
      log.printf("  using a half neighbor list\n");
    }
  } else {
    requestAtoms(atoms_lista);
  }

  // Define heavily used expressions
  sqrt2piSigma = std::sqrt(2*pi)*sigma;
  sigmaSqr2 = 2.*sigma*sigma;
  sigmaSqr = sigma*sigma;
  deltar=maxr/(nhist-1.);
  if(deltar>sigma) error("Bin size too large! Increase NHIST");
  deltaBin = std::floor(3*sigma/deltar); // 3*sigma is hard coded
  vectorX.resize(nhist);
  vectorX2.resize(nhist);
  for(unsigned i=0;i<nhist;i++){
    vectorX[i]=deltar*i;
    vectorX2[i]=vectorX[i]*vectorX[i];
  }
}

PairEntropy::~PairEntropy(){
  if (doneigh) {
     nl->printStats();
     delete nl;
  }
  if (doOutputGofr) gofrOfile.close();
  if (doOutputIntegrand) integrandOfile.close();
}

void PairEntropy::prepare(){
  if(doneigh && nl->getStride()>0){
    if(firsttime) {
      invalidateList=true;
      firsttime=false;
    } else if ( (nl->getStride()>=0) &&  (getStep()%nl->getStride()==0) ){
      invalidateList=true;
    } else if ( (nl->getStride()<0) && !(nl->isListStillGood(getPositions())) ){
      invalidateList=true;
    } else {
      invalidateList=false;
    }
  }
}

// calculator
void PairEntropy::calculate()
{
  // Define intermediate quantities
  vector<double> gofr(nhist);
  Matrix<Vector> gofrPrime(getNumberOfAtoms(),nhist);
  vector<Tensor> gofrVirial(nhist);
  // Setup parallelization
  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  if(serial){
    stride=1;
    rank=0;
  }else{
    stride=comm.Get_size();
    rank=comm.Get_rank();
  }
  // Normalization constant
  double volume=getBox().determinant();
  double density;
  if (density_given>0) density=density_given;
  else density=getNumberOfAtoms()/volume;
  double TwoPiDensity = 2*pi*density;
  double normConstantBase = TwoPiDensity*getNumberOfAtoms(); // Normalization of g(r)
  normConstantBase *= sqrt2piSigma; // Normalization of gaussian
  double invNormConstantBase = 1./normConstantBase; 
  // Calculation of g(r)
  if (doneigh && !doLowComm) {
    if(invalidateList){
      nl->update(getPositions());
    }
    // Loop over all atoms
    for(unsigned int i=0;i<nl->getNumberOfLocalAtoms();i++) {
       std::vector<unsigned> neighbors;
       unsigned index=nl->getIndexOfLocalAtom(i);
       neighbors=nl->getNeighbors(index);
       Vector position_index=getPosition(index);
       // Loop over neighbors
       for(unsigned int j=0;j<neighbors.size();j++) {  
         unsigned neighbor=neighbors[j];
         Vector distance=pbcDistance(position_index,getPosition(neighbor));
         double d2;
         if ( (d2=distance[0]*distance[0])<rcut2 && (d2+=distance[1]*distance[1])<rcut2 && (d2+=distance[2]*distance[2])<rcut2) {
           double distanceModulo=std::sqrt(d2);
           Vector distance_versor = distance / distanceModulo;
           unsigned bin=std::floor(distanceModulo/deltar);
           int minBin, maxBin; // These cannot be unsigned
           // Only consider contributions to g(r) of atoms less than n*sigma bins apart from the actual distance
           minBin=bin - deltaBin;
           if (minBin < 0) minBin=0;
           if (minBin > (nhist-1)) minBin=nhist-1;
           maxBin=bin +  deltaBin;
           if (maxBin > (nhist-1)) maxBin=nhist-1;
           for(int k=minBin;k<maxBin+1;k++) {
             double invNormKernel=invNormConstantBase/vectorX2[k];
             double dfunc;
             gofr[k] += kernel(vectorX[k]-distanceModulo,invNormKernel,dfunc);
             if (!doNotCalculateDerivatives()) {
                Vector value = dfunc * distance_versor;
                gofrPrime[index][k] += value;
                gofrPrime[neighbor][k] -= value;
                Tensor vv(value, distance);
                gofrVirial[k] += vv;
             }
           }
         }
       }
    }
  } else if (doneigh && doLowComm) {
    if(invalidateList){
      nl->update(getPositions());
    }
    // Loop over all atoms
    for(unsigned int i=0;i<nl->getNumberOfLocalAtoms();i++) {
       std::vector<unsigned> neighbors;
       unsigned index=nl->getIndexOfLocalAtom(i);
       neighbors=nl->getNeighbors(index);
       Vector position_index=getPosition(index);
       // Loop over neighbors
       for(unsigned int j=0;j<neighbors.size();j++) {  
         unsigned neighbor=neighbors[j];
         if(getAbsoluteIndex(index)==getAbsoluteIndex(neighbor)) continue;
         Vector distance=pbcDistance(position_index,getPosition(neighbor));
         double d2;
         if ( (d2=distance[0]*distance[0])<rcut2 && (d2+=distance[1]*distance[1])<rcut2 && (d2+=distance[2]*distance[2])<rcut2) {
           double distanceModulo=std::sqrt(d2);
           Vector distance_versor = distance / distanceModulo;
           unsigned bin=std::floor(distanceModulo/deltar);
           int minBin, maxBin; // These cannot be unsigned
           // Only consider contributions to g(r) of atoms less than n*sigma bins apart from the actual distance
           minBin=bin - deltaBin;
           if (minBin < 0) minBin=0;
           if (minBin > (nhist-1)) minBin=nhist-1;
           maxBin=bin +  deltaBin;
           if (maxBin > (nhist-1)) maxBin=nhist-1;
           for(int k=minBin;k<maxBin+1;k++) {
             double invNormKernel=invNormConstantBase/vectorX2[k];
             double dfunc;
             gofr[k] += kernel(vectorX[k]-distanceModulo,invNormKernel,dfunc)/2.0;
             if (!doNotCalculateDerivatives()) {
                Vector value = dfunc * distance_versor;
                gofrPrime[index][k] += value;
                Tensor vv(value/2.0, distance);
                gofrVirial[k] += vv;
             }
           }
         }
       }
    }
  } else {
    for(unsigned int i=rank;i<(getNumberOfAtoms()-1);i+=stride) {
      for(unsigned int j=i+1;j<getNumberOfAtoms();j++) {
         double d2;
         Vector distance=pbcDistance(getPosition(i),getPosition(j));
         if ( (d2=distance[0]*distance[0])<rcut2 && (d2+=distance[1]*distance[1])<rcut2 && (d2+=distance[2]*distance[2])<rcut2) {
           double distanceModulo=std::sqrt(d2);
           Vector distance_versor = distance / distanceModulo;
           unsigned bin=std::floor(distanceModulo/deltar);
           int minBin, maxBin; // These cannot be unsigned
           // Only consider contributions to g(r) of atoms less than n*sigma bins apart from the actual distance
           minBin=bin - deltaBin;
           if (minBin < 0) minBin=0;
           if (minBin > (nhist-1)) minBin=nhist-1;
           maxBin=bin +  deltaBin;
           if (maxBin > (nhist-1)) maxBin=nhist-1;
           for(int k=minBin;k<maxBin+1;k++) {
             double invNormKernel=invNormConstantBase/vectorX2[k];
             double dfunc;
             gofr[k] += kernel(vectorX[k]-distanceModulo,invNormKernel,dfunc);
             if (!doNotCalculateDerivatives()) {
                Vector value = dfunc * distance_versor;
                gofrPrime[i][k] += value;
                gofrPrime[j][k] -= value;
                Tensor vv(value, distance);
                gofrVirial[k] += vv;
             }
           }
         }
      }
    }
  }
  if(!serial){
    comm.Sum(gofr);
    if (!doNotCalculateDerivatives()) {
       if (!doLowComm) {
          comm.Sum(gofrPrime);
       }
       comm.Sum(gofrVirial);
    }
  }
  // Average g(r)
  if (doAverageGofr) {
     if (!doNotCalculateDerivatives()) error("Cannot use the AVERAGE_GOFR keyword when biasing");
     for(unsigned i=0;i<nhist;i++){
        avgGofr[i] += (gofr[i]-avgGofr[i])/( (double) iteration);
        gofr[i] = avgGofr[i];
     }
     iteration++;
  }
  // Average g(r)
  if (noPrefactor and !doNotCalculateDerivatives()) error("Cannot use the NO_PREFACTOR keyword when biasing");
  // Output of gofr
  if (doOutputGofr && (getStep()%outputStride==0)) outputGofr(gofr);
  // Find where g(r) is different from zero
  unsigned j=0;
  unsigned nhist_min=0;
  while (gofr[j]<1.e-10) {
     nhist_min=j;
     ++j;
  }
  // Construct integrand
  vector<double> logGofrX2(nhist);
  vector<double> integrand(nhist);
  for(unsigned j=0;j<nhist;j++){
    if (doReferenceGofr) {
       if (referenceGofr[j]<1.e-10) {
          // Not sure about this choice
          logGofrX2[j] = 0.;
       } else {
          logGofrX2[j] = std::log(gofr[j]/referenceGofr[j])*vectorX2[j];
       }
       if (gofr[j]<1.e-10) {
          integrand[j] = referenceGofr[j]*vectorX2[j];
       } else {
          integrand[j] = (gofr[j]*logGofrX2[j])+(-gofr[j]+referenceGofr[j])*vectorX2[j];
       }
    } else {
       logGofrX2[j] = std::log(gofr[j])*vectorX2[j];
       if (gofr[j]<1.e-10) {
          integrand[j] = vectorX2[j];
       } else {
          integrand[j] = (gofr[j]*logGofrX2[j])+(-gofr[j]+1)*vectorX2[j];
       }
    }
  }
  // Output of integrands
  if (doOutputIntegrand && (getStep()%outputStride==0)) outputIntegrand(integrand);
  // Integrate to obtain pair entropy;
  double pairEntropy;
  if (!noPrefactor) pairEntropy = -TwoPiDensity*integrate(integrand,deltar);
  else pairEntropy = integrate(integrand,deltar);
  // Construct integrand and integrate derivatives
  vector<Vector> deriv(getNumberOfAtoms());
  Tensor virial;
  if (!doNotCalculateDerivatives() ) {
    if (!doLowComm) {
       // Processors have already shared the gofrPrime
       for(unsigned int j=rank;j<getNumberOfAtoms();j+=stride) {
          vector<Vector> integrandDerivatives(nhist);
          for(unsigned k=nhist_min;k<nhist;k++){
            if (gofr[k]>1.e-10) {
              integrandDerivatives[k] = gofrPrime[j][k]*logGofrX2[k];
            }
          }
          // Integrate
          deriv[j] = -TwoPiDensity*integrate(integrandDerivatives,deltar);
       }
    } else {
       // Each processor handles only its own atoms
       for(unsigned int j=0;j<nl->getNumberOfLocalAtoms();j++) {
          unsigned index=nl->getIndexOfLocalAtom(j);
          vector<Vector> integrandDerivatives(nhist);
          for(unsigned k=nhist_min;k<nhist;k++){
            if (gofr[k]>1.e-10) {
              integrandDerivatives[k] = gofrPrime[index][k]*logGofrX2[k];
            }
          }
          // Integrate
          deriv[index] = -TwoPiDensity*integrate(integrandDerivatives,deltar);
       }
    }
    if(!serial){
      comm.Sum(deriv);
    }
    // Virial of positions
    // Construct virial integrand
    vector<Tensor> integrandVirial(nhist);
    for(unsigned j=nhist_min;j<nhist;j++){
      if (gofr[j]>1.e-10) {
        integrandVirial[j] = gofrVirial[j]*logGofrX2[j];
      }
    }
    // Integrate virial
    virial = -TwoPiDensity*integrate(integrandVirial,deltar);
    // Virial of volume
    if (density_given<0) {
      // Construct virial integrand
      vector<double> integrandVirialVolume(nhist);
      for(unsigned j=0;j<nhist;j++) {
        if (doReferenceGofr) {
           integrandVirialVolume[j] = (-gofr[j]+referenceGofr[j])*vectorX2[j];
        } else {
           integrandVirialVolume[j] = (-gofr[j]+1)*vectorX2[j];
        }
      }
      // Integrate virial
      virial += -TwoPiDensity*integrate(integrandVirialVolume,deltar)*Tensor::identity();
      }
  }
  // Assign output quantities
  for(unsigned i=0;i<deriv.size();i++) setAtomsDerivatives(i,deriv[i]);
  setValue           (pairEntropy);
  setBoxDerivatives  (virial);
}

double PairEntropy::kernel(double distance,double invNormKernel, double&der)const{
  // Gaussian function and derivative
  double result = invNormKernel*std::exp(-distance*distance/sigmaSqr2) ;
  der = -distance*result/sigmaSqr;
  return result;
}

double PairEntropy::integrate(vector<double> integrand, double delta)const{
  // Trapezoid rule
  double result = 0.;
  for(unsigned i=1;i<(integrand.size()-1);i++){
    result += integrand[i];
  }
  result += 0.5*integrand[0];
  result += 0.5*integrand[integrand.size()-1];
  result *= delta;
  return result;
}

Vector PairEntropy::integrate(vector<Vector> integrand, double delta)const{
  // Trapezoid rule
  Vector result;
  for(unsigned i=1;i<(integrand.size()-1);i++){
      result += integrand[i];
  }
  result += 0.5*integrand[0];
  result += 0.5*integrand[integrand.size()-1];
  result *= delta;
  return result;
}

Tensor PairEntropy::integrate(vector<Tensor> integrand, double delta)const{
  // Trapezoid rule
  Tensor result;
  for(unsigned i=1;i<(integrand.size()-1);i++){
      result += integrand[i];
  }
  result += 0.5*integrand[0];
  result += 0.5*integrand[integrand.size()-1];
  result *= delta;
  return result;
}

void PairEntropy::outputGofr(vector<double> gofr) {
  for(unsigned i=0;i<gofr.size();i++){
     gofrOfile.printField("r",vectorX[i]).printField("gofr",gofr[i]).printField();
  }
  gofrOfile.printf("\n");
  gofrOfile.printf("\n");
}

void PairEntropy::outputIntegrand(vector<double> integrand) {
  for(unsigned i=0;i<integrand.size();i++){
     integrandOfile.printField("r",vectorX[i]).printField("integrand",integrand[i]).printField();
  }
  integrandOfile.printf("\n");
  integrandOfile.printf("\n");
}


}
}

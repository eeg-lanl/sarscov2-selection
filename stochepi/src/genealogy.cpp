#include "genealogy.hpp"

StatePath makeStatePath(const ParticlePath & ppath) {
  StatePath spath;
  for ( auto & particle : ppath ) {
    spath.splice(spath.end(), particle.getPath());
  }
  return spath;
}




void Genealogy::update(const std::vector<Particle*> & particle_vec,
    const std::map<int, int> & counts) {
  // add sampled particles to the ancestors
  ParticleIdxMap pim; AncestorIdxMap aim;
  // k keeps track of the indices of the new particles
  int k = 0;
  for ( auto & [idx, count] : counts ) {
    // add Particle to the Particle map, but only if count > 0
    if ( count > 0 ) {
      pim.emplace(idx, *(particle_vec[idx]));
    }
    /* the new particles are created by walking through the sorted indices
     * and adding c copies to the new list.
     * Because a map is sorted, we can add c copies of idx as values to the
     * AncestorIdxMap
     */
    for ( int c = 0; c < count; ++c ) {
      aim[k] = idx; ++k;
    }
  }
  // add pim and aim to the list (using std::move)
  particle_idx_maps.push_back(std::move(pim));
  ancestries.push_back(std::move(aim));
  // walk through ancestors, and remove extinct lineages
  auto prit = particle_idx_maps.rbegin();
  auto arit = ancestries.rbegin(); arit++;
  while ( prit != particle_idx_maps.rend() && arit != ancestries.rend() ) {
    auto & pim = *prit; auto & aim = *arit;
    // remove keys from ancestry that are not in particle map
    auto ait = aim.begin();
    while ( ait != aim.end() ) {
      // try to find it->first in keys if pim
      if ( pim.find(ait->first) == pim.end() ) {
        // remove it from ancestry
        ait = aim.erase(ait);
      } else {
        // just increase the iterator
        ++ait;
      }
    }
    // now remove particles that are not pointed at in the ancestry
    ++prit;
    if ( prit == particle_idx_maps.rend() ) {
      throw std::logic_error("particle and ancestry lists have unequal size" + RIGHT_HERE);
    } // else the iterator is valid
    auto & ppim = *prit;
    auto pit = ppim.begin();
    while ( pit != ppim.end() ) {
      // find pit->first in values of aim
      int idx = pit->first;
      auto pred = [idx](auto & kv){return kv.second == idx;};
      if ( std::find_if(aim.begin(), aim.end(), pred) == aim.end() ) {
        // remove particle
        pit = ppim.erase(pit);
      } else {
        // just increase the iterator
        ++pit;
      }
    }
    if ( ppim.empty() ) {
      throw std::logic_error("Particle Index map should never become empty" + RIGHT_HERE);
    }
    // finally, increase the ancestry iterator
    ++arit;
  }
}

bool Genealogy::empty() const {
  return ( particle_idx_maps.empty() || ancestries.empty() );
}

void Genealogy::clear() {
  particle_idx_maps.clear();
  ancestries.clear();
}

ParticlePath Genealogy::sampleParticlePath(Rng & rng) const {
  ParticlePath path;
  if ( empty() ) {
    return path;
  } // else: the Genealogy is non-empty
  // FORWARD SAMPLING
  auto pit = particle_idx_maps.begin();
  auto ait = ancestries.begin();
  // uniformly sample a random Particle from pim
  auto & pim = *pit;
  int r = rng.integer(pim.size());
  auto it = pim.begin(); std::advance(it, r);
  int idx = it->first;
  // now loop through the lists
  while ( pit != particle_idx_maps.end() && ait != ancestries.end() ) {
    // get aliases
    auto & pim = *pit; auto & aim = *ait;
    // add the particle to the path
    path.push_back(pim.at(idx));
    // sample a descendent of the Particle
    auto cpred = [idx](auto & x){return x.second == idx;};
    // but first count the number of descendents
    int num_desc = std::count_if(aim.begin(), aim.end(), cpred);
    if ( num_desc == 0 ) {
      throw std::logic_error("Particles in the Genealogy should have descendents" + RIGHT_HERE);
    }
    int r = rng.integer(num_desc);
    auto fpred = [&r, idx](auto & x){
      if ( x.second == idx ) r--;
      return ( r < 0 );
    };
    auto it = std::find_if(aim.begin(), aim.end(), fpred);
    idx = it->first;
    // and increment the iterators
    ++pit; ++ait;
  }
  return path;
}

StatePath Genealogy::sampleStatePath(Rng & rng) const {
  return makeStatePath(sampleParticlePath(rng));
}


void Genealogy::print(std::ostream & os) const {
  auto ait = ancestries.begin();
  auto pit = particle_idx_maps.begin();
  while ( ait != ancestries.end() && pit != particle_idx_maps.end() ) {
    std::cout << "particles:" << std::endl;
    auto & pim = *pit;
    for ( auto & [idx, p] : pim ) {
      (void) p;
      std::cout << idx << ", ";
    }
    std::cout << "\nancestries:" << std::endl;
    auto & aim = *ait;
    for ( auto & [idx1, idx2] : aim ) {
      std::cout << idx1 << "->" << idx2 << ", ";
    }
    std::cout << std::endl;
    // increment iterators...
    ++ait; ++pit;
  }
}
